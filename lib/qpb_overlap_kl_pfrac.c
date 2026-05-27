#include <string.h>
#include <qpb_types.h>
#include <qpb_errors.h>
#include <qpb_globals.h>
#include <qpb_alloc.h>
#include <qpb_spinor_field.h>
#include <qpb_spinor_linalg.h>
#include <qpb_gauge_field.h>
#include <qpb_comm_halo_spinor_field.h>
#include <qpb_comm_halo_gauge_field.h>
#include <qpb_timebc_set_gauge_field.h>
#include <qpb_dslash_wrappers.h>
#include <qpb_stop_watch.h>
#include <qpb_kl_defs.h>
#include <qpb_mscongrad.h>
#include <math.h>

/*
  ov_temp_vecs layout (16 vectors total):

   Apply scratch (transient; clobbered on every call to the apply functions):
    [ 0]  sum vector for the partial-fraction accumulation
    [ 1]  z vector  (D·sum in apply_overlap)
    [ 2]  γ5 scratch in apply_conj_overlap

   Outer BiCGStab persistent state:
    [ 3]  r0_hat                         [ 7]  v = D_ov^outer · z_pc
    [ 4]  r                              [ 8]  y_pc  = K^{-1} · p
    [ 5]  p                              [ 9]  z_pc  = K^{-1} · s
    [ 6]  u = D_ov^outer · y_pc

   Preconditioner-solver scratch (one set of slots, two possible owners):
     n_outer ≥ 2 → inner BiCGStab on D_ov^(n_outer-1):
       [10]  r0_hat_in   [11] r_in   [12] p_in   [13] u_in   [14] v_in   [15] (spare)
     n_outer  = 1 → legacy CGNE on the shifted kernel (kept verbatim):
       [10..15]  r, p, w, y, bprime, bpp  (same six slots as legacy nPrec=0)

   The two preconditioner paths never coexist at runtime — n_outer is fixed
   at init — so they safely share these six slots.
*/

#define OVERLAP_NUMB_TEMP_VECS 16
#define MSCG_NUMB_TEMP_VECS    20

typedef enum {
  LEVEL_OUTER = 0,
  LEVEL_PREC  = 1,
  LEVEL_COUNT = 2
} overlap_level_t;

static qpb_spinor_field ov_temp_vecs  [OVERLAP_NUMB_TEMP_VECS];
static qpb_spinor_field mscg_temp_vecs[MSCG_NUMB_TEMP_VECS];

static qpb_overlap_params ov_params;

/* Per-level partial-fraction data.
   For n_outer = 1, LEVEL_PREC is unused (preconditioner is the shifted
   kernel, not an overlap). */
static int          kl_order      [LEVEL_COUNT];
static qpb_double  *shifts        [LEVEL_COUNT];
static qpb_double  *numerators    [LEVEL_COUNT];
static qpb_double   constant_term [LEVEL_COUNT];

/* MSCG control: per-level tolerance, shared iteration cap. */
static qpb_double MS_solver_precision[LEVEL_COUNT];
static int        MS_maximum_solver_iterations;

/* Preconditioner-solver control (one set of knobs; semantics adapt to path):
     n_outer = 1 → drives legacy CGNE-on-shifted-kernel
     n_outer ≥ 2 → drives the new inner BiCGStab on D_ov^(n_prec) */
static qpb_double prec_solver_epsilon;
static int        prec_solver_max_iter;


void
qpb_overlap_kl_pfrac_init(void * gauge, qpb_clover_term clover,
          enum qpb_kl_classes kl_class, int kl_iters, qpb_double rho,
          qpb_double c_sw, qpb_double mass, qpb_double scaling_factor,
          qpb_double ms_epsilon, qpb_double prec_ms_epsilon, int ms_max_iters,
          qpb_double prec_epsilon, int prec_max_iter)
{
  if(ov_params.initialized == QPB_OVERLAP_INITIALIZED)
    return;

  qpb_comm_halo_spinor_field_init();
  for(int i=0; i<OVERLAP_NUMB_TEMP_VECS; i++) {
    ov_temp_vecs[i] = qpb_spinor_field_init();
    qpb_spinor_field_set_zero(ov_temp_vecs[i]);
  }
  for(int i=0; i<MSCG_NUMB_TEMP_VECS; i++) {
    mscg_temp_vecs[i] = qpb_spinor_field_init();
    qpb_spinor_field_set_zero(mscg_temp_vecs[i]);
  }


  qpb_gauge_field gauge_bc;
  if(which_dslash_op == QPB_DSLASH_STANDARD)
  {
    gauge_bc = qpb_gauge_field_init();
    qpb_timebc_set_gauge_field(gauge_bc, *(qpb_gauge_field *)gauge,\
                  problem_params.timebc);
    ov_params.gauge_ptr = qpb_alloc(sizeof(qpb_gauge_field));
    memcpy(ov_params.gauge_ptr, &gauge_bc, sizeof(qpb_gauge_field));
  }
  else
  {
    ov_params.gauge_ptr = gauge;
  }

  ov_params.c_sw = c_sw;
  ov_params.rho = rho;
  ov_params.m_bare = -rho; // Kernel operator bare mass
  ov_params.mass = mass;
  ov_params.clover = clover;

  switch(which_dslash_op)
  {
  case QPB_DSLASH_BRILLOUIN:
    if(c_sw)
    {
      ov_params.g5_dslash_op = &qpb_gamma5_clover_bri_dslash;
      ov_params.dslash_op = &qpb_clover_bri_dslash;
    }
    else
    {
      ov_params.g5_dslash_op = &qpb_gamma5_bri_dslash;
      ov_params.dslash_op = &qpb_bri_dslash;
    }
    break;
  case QPB_DSLASH_STANDARD:
    if(c_sw)
    {
      ov_params.g5_dslash_op = &qpb_gamma5_clover_dslash;
      ov_params.dslash_op = &qpb_clover_dslash;
    }
    else
    {
      ov_params.g5_dslash_op = &qpb_gamma5_dslash;
      ov_params.dslash_op = &qpb_dslash;
    }
    break;
  }
  ov_params.initialized = QPB_OVERLAP_INITIALIZED;

  if(kl_iters < 1) {
    error("qpb_overlap_kl_pfrac_init: kl_iters must be ≥ 1, got %d\n",
          kl_iters);
    exit(QPB_PARAMETERS_ERROR);
  }

  /* Levels: outer = kl_iters, prec = kl_iters - 1. Prec level is only
     meaningfully populated when n_outer ≥ 2; at n_outer = 1 the
     preconditioner falls back to the shifted kernel. */
  kl_order[LEVEL_OUTER] = kl_iters;
  kl_order[LEVEL_PREC]  = kl_iters - 1;

  MS_solver_precision[LEVEL_OUTER] = ms_epsilon;
  MS_solver_precision[LEVEL_PREC]  = prec_ms_epsilon;
  MS_maximum_solver_iterations     = ms_max_iters;

  prec_solver_epsilon  = prec_epsilon;
  prec_solver_max_iter = prec_max_iter;

  /* Populate partial-fraction tables per level. Skip LEVEL_PREC when its
     order is 0 — the kernel path doesn't use shifts/numerators. */
  for(int lvl = 0; lvl < LEVEL_COUNT; lvl++) {
    int n = kl_order[lvl];
    if(n == 0) { shifts[lvl] = NULL; numerators[lvl] = NULL;
                 constant_term[lvl] = 0.; continue; }

    shifts    [lvl] = qpb_alloc(sizeof(qpb_double) * n);
    numerators[lvl] = qpb_alloc(sizeof(qpb_double) * n);
    constant_term[lvl] = 1.0 / ((qpb_double)(2*n + 1));

    for(int i = 0; i < n; i++) {
      qpb_double t = M_PI * (i + 0.5) * constant_term[lvl];
      shifts    [lvl][i] = pow(tan(t), 2);
      numerators[lvl][i] = 2 * constant_term[lvl] / powl(cos(t), 2);
    }

    if(scaling_factor != 1.0) {
      constant_term[lvl] *= 1.0 / sqrt(scaling_factor);
      for(int i = 0; i < n; i++) {
        numerators[lvl][i] *= sqrt(scaling_factor);
        shifts    [lvl][i] *= scaling_factor;
      }
    }
  }

  /* MSCG workspace sized for the larger of the two orders = n_outer. */
  qpb_mscongrad_init(kl_order[LEVEL_OUTER]);
}


void
qpb_overlap_kl_pfrac_finalize(void)
{
  qpb_comm_halo_spinor_field_finalize();
  for(int i=0; i<OVERLAP_NUMB_TEMP_VECS; i++)
    qpb_spinor_field_finalize(ov_temp_vecs[i]);
  for(int i=0; i<MSCG_NUMB_TEMP_VECS; i++)
    qpb_spinor_field_finalize(mscg_temp_vecs[i]);

  if(which_dslash_op == QPB_DSLASH_STANDARD)
    qpb_gauge_field_finalize(*(qpb_gauge_field *)ov_params.gauge_ptr);

  ov_params.initialized = 0;
  qpb_mscongrad_finalize(kl_order[LEVEL_OUTER]);

  for(int lvl = 0; lvl < LEVEL_COUNT; lvl++) {
    if(shifts    [lvl]) free(shifts    [lvl]);
    if(numerators[lvl]) free(numerators[lvl]);
  }
}



INLINE void
D_op(qpb_spinor_field y, qpb_spinor_field x)
{
  /* Implements D - rho*I */

  void *dslash_args[4];

  dslash_args[0] = ov_params.gauge_ptr;
  dslash_args[1] = &ov_params.m_bare;
  dslash_args[2] = &ov_params.clover;
  dslash_args[3] = &ov_params.c_sw;

  ov_params.dslash_op(y, x, dslash_args);

  return;
}


static void
apply_gamma5_sign(overlap_level_t lvl,
                  qpb_spinor_field y, qpb_spinor_field x)
{
  
  /* γ5·sign(X)·x at the partial-fraction order of this level. */
  qpb_spinor_field sum = ov_temp_vecs[0];
  
  int n = kl_order[lvl];
  if(n == 0) { error("apply_gamma5_sign: level %d has KL order 0\n", lvl);
               exit(QPB_PARAMETERS_ERROR); }

  qpb_spinor_field yMS[n];
  for(int s = 0; s < n; s++) {
    yMS[s] = mscg_temp_vecs[s];
    qpb_spinor_field_set_zero(yMS[s]);
  }

  qpb_double kernel_mass  = ov_params.m_bare;
  qpb_double kernel_kappa = 1.0 / (2*kernel_mass + 8.0);

  qpb_mscongrad(yMS, x, ov_params.gauge_ptr, ov_params.clover, kernel_kappa,
                n, shifts[lvl], ov_params.c_sw,
                MS_solver_precision[lvl], MS_maximum_solver_iterations);

  qpb_spinor_ax(sum, (qpb_complex){constant_term[lvl], 0.}, x);
  for(int s = 0; s < n; s++)
    qpb_spinor_axpy(sum, (qpb_complex){numerators[lvl][s], 0.},
                    yMS[s], sum);

  D_op(y, sum);
}


static void
apply_overlap(overlap_level_t lvl,
              qpb_spinor_field y, qpb_spinor_field x)
{
  /* D_ov,m^(lvl)·x = (ρ + a m/2)·x + (ρ - a m/2)·γ5·sign(X)·x .
     Preconditioner mass = overlap mass (no separate parameter). */
  qpb_spinor_field z = ov_temp_vecs[1];

  qpb_double m   = ov_params.mass;
  qpb_double rho = ov_params.rho;
  qpb_complex a  = { rho + 0.5*m, 0. };
  qpb_complex b  = { rho - 0.5*m, 0. };

  apply_gamma5_sign(lvl, z, x);
  qpb_spinor_axpby(y, a, x, b, z);
}

static void
apply_conj_overlap(overlap_level_t lvl,
                   qpb_spinor_field y, qpb_spinor_field x)
{
  /* D_ov† = γ5 · D_ov · γ5 (γ5-hermiticity). */
  qpb_spinor_field w = ov_temp_vecs[2];

  qpb_spinor_gamma5(w, x);
  apply_overlap(lvl, y, w);
  qpb_spinor_gamma5(y, y);
}

/* Public single-arg wrappers (callers from main programs). */
void qpb_overlap_kl_pfrac          (qpb_spinor_field y, qpb_spinor_field x)
                                   { apply_overlap     (LEVEL_OUTER, y, x); }
void qpb_conjugate_overlap_kl_pfrac(qpb_spinor_field y, qpb_spinor_field x)
                                   { apply_conj_overlap(LEVEL_OUTER, y, x); }


static int
preconditioner_bicgstab(qpb_spinor_field x, qpb_spinor_field b)
{
  qpb_spinor_field r0_hat = ov_temp_vecs[10];
  qpb_spinor_field r      = ov_temp_vecs[11];
  qpb_spinor_field p      = ov_temp_vecs[12];
  qpb_spinor_field u      = ov_temp_vecs[13];
  qpb_spinor_field v      = ov_temp_vecs[14];

  qpb_double res_norm, b_norm;
  qpb_complex_double alpha = {1, 0}, omega = {1, 0};
  qpb_complex_double beta, gamma, rho, zeta;
  int iters;

  qpb_spinor_xdotx(&b_norm, b);

  qpb_spinor_field_set_zero(x);
  qpb_spinor_xeqy(r, b);            /* r0 = b - D_ov^prec·x0 = b */
  qpb_spinor_xeqy(r0_hat, r);

  qpb_spinor_xdotx(&res_norm, r);
  gamma.re = res_norm; gamma.im = 0.;
  rho = gamma;

  for(iters = 1; iters < prec_solver_max_iter; iters++) {
    if(res_norm / b_norm <= prec_solver_epsilon * prec_solver_epsilon)
      break;

    qpb_spinor_xdoty(&gamma, r0_hat, r);
    beta = CMUL(CDEV(gamma, rho), CDEV(alpha, omega));

    omega = CNEGATE(omega);
    qpb_spinor_axpy(p, omega, u, p);     /* p -= omega·u                 */
    qpb_spinor_axpy(p, beta, p, r);      /* p  = beta·p + r              */

    apply_overlap(LEVEL_PREC, u, p);     /* u  = D_ov^prec · p           */

    qpb_spinor_xdoty(&beta, r0_hat, u);
    rho   = gamma;
    alpha = CDEV(rho, beta);

    alpha = CNEGATE(alpha);
    qpb_spinor_axpy(r, alpha, u, r);     /* r -= alpha·u   (r ≡ s now)   */

    apply_overlap(LEVEL_PREC, v, r);     /* v  = D_ov^prec · s           */

    qpb_spinor_xdoty(&zeta, v, r);
    qpb_spinor_xdotx(&beta.re, v); beta.im = 0;
    omega = CDEV(zeta, beta);

    alpha = CNEGATE(alpha);
    qpb_spinor_axpy(x, alpha, p, x);     /* x += alpha·p                 */
    qpb_spinor_axpy(x, omega, r, x);     /* x += omega·s   (s is in r)   */

    omega = CNEGATE(omega);
    qpb_spinor_axpy(r, omega, v, r);     /* r -= omega·v                 */
    omega = CNEGATE(omega);

    qpb_spinor_xdotx(&res_norm, r);
  }

  return iters;            /* caller can compare against prec_solver_max_iter */
}

// Keep verbatim from the legacy file:

static void
M_kernel_op(qpb_spinor_field y, qpb_spinor_field x)
{
  /* Implements M = c·I + r·D ,
     with c = ρ + a*m/2 (preconditioner center)
     and  r = ρ - a*m/2 (preconditioner radius).
     This is the simple kernel-based preconditioner used for nPrec=0. */

  qpb_spinor_field w = ov_temp_vecs[0];

  qpb_double rho = ov_params.rho;
  qpb_double preconditioner_mass = ov_params.mass;

  qpb_complex preconditioner_center = {rho + 0.5*preconditioner_mass, 0.};
  qpb_complex preconditioner_radius = {rho - 0.5*preconditioner_mass, 0.};

  D_op(w, x);
  qpb_spinor_axpby(y, preconditioner_center, x, preconditioner_radius, w);
  
  return;
}


static void
M_kernel_conj_op(qpb_spinor_field y, qpb_spinor_field x)
{
  /* Implements M† = γ5·M·γ5 ,
     valid because of the γ5-hermiticity of D: γ5·D·γ5 = D†. */

  qpb_spinor_field z = ov_temp_vecs[1];

  qpb_spinor_gamma5(y, x);
  M_kernel_op(z, y);
  qpb_spinor_gamma5(y, z);
  
  return;
}


static int
kernel_preconditioner_CGNE(qpb_spinor_field x, qpb_spinor_field b)
{
  /* CGNE on the shifted kernel M = c·I + r·D; used by
  qpb_bicgstab_overlap_kl_pfrac only when n_outer == 1. */

  qpb_spinor_field r      = ov_temp_vecs[10];
  qpb_spinor_field p      = ov_temp_vecs[11];
  qpb_spinor_field w      = ov_temp_vecs[12];
  qpb_spinor_field y      = ov_temp_vecs[13];
  qpb_spinor_field bprime = ov_temp_vecs[14];
  qpb_spinor_field bpp    = ov_temp_vecs[15];

  int iters = 0;

  /* All scalars are real -- K†K is HPD */
  qpb_double bpp_norm, res_norm, new_res_norm, omega, alpha, beta;
  
  qpb_spinor_xeqy(bprime, b);

  /* bpp = K†·bprime  (RHS of the normal-equations system K†K·x = K†·bprime) */
  M_kernel_conj_op(bpp, bprime);
  qpb_spinor_xdotx(&bpp_norm, bpp);

  /* --- x0 = 0 --- */
  qpb_spinor_field_set_zero(x);
  qpb_spinor_xeqy(r, bpp);

  /* gamma_0 = ||r0||^2 */
  qpb_spinor_xdotx(&res_norm, r);

  /* p0 = r0 */
  qpb_spinor_xeqy(p, r);

  for(iters = 1; iters < prec_solver_max_iter; iters++)
  {
    /* Stopping on relative residual of the normal-equations system:
       ||r||^2 / ||K†·bprime||^2 <= prec_solver_epsilon^2 */
    if(res_norm / bpp_norm <= prec_solver_epsilon * prec_solver_epsilon)
      break;

    /* Apply K†K to p: w = K·p,  y = K†·w = K†K·p */
    M_kernel_op(w, p);
    M_kernel_conj_op(y, w);

    /* omega = p†·K†K·p = ||K·p||^2 = ||w||^2  (real, positive) */
    qpb_spinor_xdotx(&omega, w);

    /* alpha = ||r||^2 / (p†·K†K·p)  (real, positive) */
    alpha = res_norm / omega;

    /* x += alpha·p */
    qpb_spinor_axpy(x, (qpb_complex) {alpha, 0.}, p, x);

    /* r -= alpha·K†K·p = alpha·y */
    qpb_spinor_axpy(r, (qpb_complex) {-alpha, 0.}, y, r);
    
    /* new_res_norm = ||r_{k+1}||^2 */
    qpb_spinor_xdotx(&new_res_norm, r);

    /* beta = ||r_{k+1}||^2 / ||r_k||^2 */
    beta = new_res_norm / res_norm;

    /* p_{k+1} = r_{k+1} + beta·p_k */
    qpb_spinor_axpy(p, (qpb_complex) {beta, 0.}, p, r);

    res_norm = new_res_norm;
  }

  return iters;
}


// M_kernel_op and M_kernel_conj_op qpb_preconditioner_CGNE — rename to
// kernel_preconditioner_CGNE for clarity and have it use slots
// [10..15]. Its body doesn't change.

int
qpb_bicgstab_overlap_kl_pfrac(qpb_spinor_field x, qpb_spinor_field b,
                              qpb_double epsilon, int max_iter)
{
  qpb_spinor_field r0_hat = ov_temp_vecs[3];
  qpb_spinor_field r      = ov_temp_vecs[4];
  qpb_spinor_field p      = ov_temp_vecs[5];
  qpb_spinor_field u      = ov_temp_vecs[6];
  qpb_spinor_field v      = ov_temp_vecs[7];
  qpb_spinor_field y_pc   = ov_temp_vecs[8];   /* K^{-1}·p */
  qpb_spinor_field z_pc   = ov_temp_vecs[9];   /* K^{-1}·s */

  int n_reeval = 100, n_echo = 100, iters = 0;
  qpb_double res_norm, b_norm;
  qpb_complex_double alpha = {1, 0}, omega = {1, 0};
  qpb_complex_double beta, gamma, rho, zeta;

  /* Select the preconditioner path once, outside the loop. */
  int n_outer = kl_order[LEVEL_OUTER];
  enum { PREC_NONE, PREC_KERNEL, PREC_BICGSTAB } prec_path;
  if(prec_solver_max_iter == 0)   prec_path = PREC_NONE;
  else if(n_outer == 1)           prec_path = PREC_KERNEL;
  else                            prec_path = PREC_BICGSTAB;

  qpb_spinor_xdotx(&b_norm, b);
  qpb_spinor_field_set_zero(x);
  qpb_spinor_xeqy(r, b);                 /* r0 = b - D_ov·x0 = b */
  qpb_spinor_xeqy(r0_hat, r);
  qpb_spinor_xdotx(&res_norm, r);
  gamma.re = res_norm; gamma.im = 0;
  rho = gamma;

  qpb_double t = qpb_stop_watch(0);
  for(iters = 1; iters < max_iter; iters++) {
    if(res_norm / b_norm <= epsilon) break;

    qpb_spinor_xdoty(&gamma, r0_hat, r);
    beta = CMUL(CDEV(gamma, rho), CDEV(alpha, omega));

    omega = CNEGATE(omega);
    qpb_spinor_axpy(p, omega, u, p);
    qpb_spinor_axpy(p, beta, p, r);

    /* y_pc = K^{-1}·p */
    switch(prec_path) {
    case PREC_NONE:     qpb_spinor_xeqy(y_pc, p);                  break;
    case PREC_KERNEL:   kernel_preconditioner_CGNE(y_pc, p);       break;
    case PREC_BICGSTAB: preconditioner_bicgstab   (y_pc, p);       break;
    }

    apply_overlap(LEVEL_OUTER, u, y_pc);          /* u = D_ov · y_pc */

    qpb_spinor_xdoty(&beta, r0_hat, u);
    rho   = gamma;
    alpha = CDEV(rho, beta);

    alpha = CNEGATE(alpha);
    qpb_spinor_axpy(r, alpha, u, r);              /* r ≡ s now */

    /* z_pc = K^{-1}·s */
    switch(prec_path) {
    case PREC_NONE:     qpb_spinor_xeqy(z_pc, r);                  break;
    case PREC_KERNEL:   kernel_preconditioner_CGNE(z_pc, r);       break;
    case PREC_BICGSTAB: preconditioner_bicgstab   (z_pc, r);       break;
    }

    apply_overlap(LEVEL_OUTER, v, z_pc);          /* v = D_ov · z_pc */

    qpb_spinor_xdoty(&zeta, v, r);
    qpb_spinor_xdotx(&beta.re, v); beta.im = 0;
    omega = CDEV(zeta, beta);

    alpha = CNEGATE(alpha);
    qpb_spinor_axpy(x, alpha, y_pc, x);           /* x += alpha · y_pc */
    qpb_spinor_axpy(x, omega, z_pc, x);           /* x += omega · z_pc */

    if(iters % n_reeval == 0) {
      apply_overlap(LEVEL_OUTER, u, x);
      qpb_spinor_xmy(r, b, u);
    } else {
      omega = CNEGATE(omega);
      qpb_spinor_axpy(r, omega, v, r);
      omega = CNEGATE(omega);
    }

    qpb_spinor_xdotx(&res_norm, r);
    if(iters % n_echo == 0)
      print(" \t iters = %8d, res = %e\n", iters, res_norm / b_norm);
  }
  t = qpb_stop_watch(t);

  /* Final explicit residual */
  apply_overlap(LEVEL_OUTER, u, x);
  qpb_spinor_xmy(r, b, u);
  qpb_spinor_xdotx(&res_norm, r);

  if(iters == max_iter) {
    error(" Preconditioned BiCGStab *did not* converge after %d iters\n", iters);
    error(" residual = %e, relative = %e, t = %g sec\n",
          res_norm, res_norm / b_norm, t);
    return -1;
  }
  print(" \tAfter %d iters, preconditioned BiCGStab converged, res = %e, "
        "relative = %e, t = %g sec\n",
        iters, res_norm, res_norm / b_norm, t);
  return iters;
}