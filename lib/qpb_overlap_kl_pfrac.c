/* ================================================================
   qpb_toy_preconditioned_cgne.c

   PURPOSE: Minimal sanity check of the preconditioned CGNE machinery.

   Outer system:     D_toy · x = b
   Outer solver:     CGNE (normal equations D†D · x = D†b)
   Preconditioner:   SAME D_toy, solved via a looser-tolerance inner CGNE

   where D_toy ≡ D_Bri + m·I  (massive Brillouin operator, mass = overlap mass).

   Because the preconditioner is IDENTICAL to the outer operator, M^{-1}A = I
   and the condition number is 1.  With a sufficiently tight inner solve the
   outer CG should converge in O(1) iterations.  Any deviation from this is
   a clear signal of a bug in the preconditioned-CG loop itself.

   γ5-Hermiticity:  D†_toy = γ5 · D_toy · γ5  (real mass m, same as D_op).

   Temp-vector slot map (toy_temp_vecs[]):
     [0]          scratch for Dconj_toy_op  (γ5·x intermediate)
     [1..6]       inner CGNE:  r, p, z, w, y, bprime
     [7..13]      outer CGNE:  r, p, z, y, w, bprime, s
   ================================================================ */

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
#include <math.h>


#define TOY_NUMB_TEMP_VECS 14

static qpb_spinor_field toy_temp_vecs[TOY_NUMB_TEMP_VECS];

// /* Parameters needed to apply D_toy = D_Bri + m·I */
// static void        *toy_gauge_ptr;
// static qpb_clover_term toy_clover;
// static qpb_double   toy_c_sw;
// static qpb_double   toy_mass;    /* the "overlap" mass m; D_toy = D_Bri + m·I */

static qpb_overlap_params ov_params;

static int KL_diagonal_order;
static qpb_double MS_solver_precision;
static int MS_maximum_solver_iterations;

static qpb_double prec_epsilon;
static int prec_max_iter;


/* ================================================================
   Init / Finalize
   ================================================================ */

void
qpb_overlap_kl_pfrac_init(void * gauge, qpb_clover_term clover, \
          enum qpb_kl_classes kl_class, int kl_iters, qpb_double rho, \
          qpb_double c_sw, qpb_double mass, qpb_double scaling_factor, \
          qpb_double ms_epsilon, int ms_max_iter)
{

  if(ov_params.initialized != QPB_OVERLAP_INITIALIZED)
  {

  // if(!toy_initialized)
  // {
    qpb_comm_halo_spinor_field_init();

    for(int i = 0; i < TOY_NUMB_TEMP_VECS; i++)
    {
      toy_temp_vecs[i] = qpb_spinor_field_init();
      qpb_spinor_field_set_zero(toy_temp_vecs[i]);
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

    // toy_gauge_ptr = gauge;
    // toy_clover    = clover;
    // toy_c_sw      = c_sw;
    // toy_mass      = mass;

    // toy_initialized = 1;

    
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

    KL_diagonal_order = kl_iters;
    MS_solver_precision = ms_epsilon;
    MS_maximum_solver_iterations = ms_max_iter;

    prec_epsilon = ms_epsilon;
    prec_max_iter = ms_max_iter;

    print(" [toy] D_toy = D_Bri + m·I,  m = %g\n", ov_params.mass);
    print(" [toy] Preconditioner = same D_toy (different tolerance)\n");
    print(" [toy] Expected outer CG iters ≈ 1 with exact inner solve\n");
  }
  return;
}


void
qpb_overlap_kl_pfrac_finalize()
{
  qpb_comm_halo_spinor_field_finalize();

  for(int i = 0; i < TOY_NUMB_TEMP_VECS; i++)
    qpb_spinor_field_finalize(toy_temp_vecs[i]);

  if(which_dslash_op == QPB_DSLASH_STANDARD)
    qpb_gauge_field_finalize(*(qpb_gauge_field *)ov_params.gauge_ptr);

  ov_params.initialized = 0;

  // toy_initialized = 0;
  return;
}


/* ================================================================
   Operator definitions
   ================================================================ */

/* D_toy_op: y = (D_Bri + m·I) x
   Calls the Brillouin dslash with bare mass = toy_mass.
   The dslash convention is (D + m_bare)·x, so m_bare = toy_mass here.
   If c_sw != 0, the clover-improved Brillouin dslash is used instead. */
INLINE void
D_toy_op(qpb_spinor_field y, qpb_spinor_field x)
{
  void *dslash_args[4];

  dslash_args[0] = ov_params.gauge_ptr;
  dslash_args[1] = &ov_params.m_bare;
  dslash_args[2] = &ov_params.clover;
  dslash_args[3] = &ov_params.c_sw;

  ov_params.dslash_op(y, x, dslash_args);

  return;
}


/* Dconj_toy_op: y = D†_toy x = γ5 · D_toy · γ5 · x
   Uses toy_temp_vecs[0] as scratch.
   Safety: slot [0] is only touched here; it is never live during
   the inner or outer CG loops. */
INLINE void
Dconj_toy_op(qpb_spinor_field y, qpb_spinor_field x)
{
  qpb_spinor_field tmp = toy_temp_vecs[0];  /* scratch for γ5·x */

  qpb_spinor_gamma5(tmp, x);     /* tmp = γ5·x            */
  D_toy_op(y, tmp);              /* y   = D_toy·(γ5·x)    */
  qpb_spinor_gamma5(y, y);       /* y   = γ5·D_toy·γ5·x   */

  return;
}


/* ================================================================
   Inner solver: unpreconditioned CGNE for  D_toy · x = b
   Solves the normal equations  D†_toy · D_toy · x = D†_toy · b.

   This is the function that acts as the preconditioner for the outer
   loop.  It is called once per outer iteration with z (the outer
   CGNE residual) as right-hand side.  The looser its tolerance, the
   cheaper it is — but making it too loose will stall the outer loop
   (see Section 9.4 of Saad on flexible variants).

   Silent: no print statements.  Returns iteration count or -1.

   Temp slots used: [1..6]
     [1] r      residual of original system  b - D_toy·x
     [2] p      search direction             (= z_0 initially)
     [3] z      normal-equations residual    D†_toy · r
     [4] w      D_toy · p
     [5] y      D†_toy · w  = (D†D) · p
     [6] bprime D†_toy · b
   ================================================================ */
int
qpb_toy_inner_CG(qpb_spinor_field x, qpb_spinor_field b,
                 qpb_double epsilon, int max_iter)
{
  qpb_spinor_field r      = toy_temp_vecs[1];
  qpb_spinor_field p      = toy_temp_vecs[2];
  qpb_spinor_field z      = toy_temp_vecs[3];
  qpb_spinor_field w      = toy_temp_vecs[4];
  qpb_spinor_field y      = toy_temp_vecs[5];
  qpb_spinor_field bprime = toy_temp_vecs[6];

  int n_reeval = 100;
  int iters    = 0;

  qpb_double res_norm, true_res_norm, b_norm, bprime_norm;
  qpb_complex_double alpha = {1, 0}, omega = {1, 0};
  qpb_complex_double beta, gamma;

  /* ||b||  – denominator for stopping criterion */
  qpb_spinor_xdotx(&b_norm, b);
  true_res_norm = b_norm;

  /* bprime = D†_toy · b */
  Dconj_toy_op(bprime, b);
  qpb_spinor_xdotx(&bprime_norm, bprime);

  /* x0 = 0 */
  qpb_spinor_field_set_zero(x);

  /* r0 = b  (exact since x0 = 0),  z0 = bprime */
  qpb_spinor_xeqy(r, b);
  qpb_spinor_xeqy(z, bprime);

  /* gamma_0 = ||z0||^2  (real, positive) */
  qpb_spinor_xdotx(&gamma.re, z);
  gamma.im = 0.;

  /* p0 = z0 */
  qpb_spinor_xeqy(p, z);

  for(iters = 1; iters < max_iter; iters++)
  {
    /* Stopping criterion on true residual of original system D_toy·x = b */
    if(true_res_norm / b_norm <= epsilon)
      break;

    /* w = D_toy · p */
    D_toy_op(w, p);
    /* y = D†_toy · w = (D†D_toy) · p */
    Dconj_toy_op(y, w);

    /* omega = p† · (D†D) · p = ||D·p||^2 = ||w||^2  (real, positive) */
    qpb_spinor_xdoty(&omega, p, y);

    /* alpha = gamma / omega */
    alpha = CDEV(gamma, omega);

    /* x += alpha · p */
    qpb_spinor_axpy(x, alpha, p, x);

    /* Update r and z */
    if(iters % n_reeval == 0)
    {
      /* Full recomputation to suppress round-off accumulation */
      D_toy_op(w, x);
      qpb_spinor_xmy(r, b, w);
      Dconj_toy_op(y, w);
      qpb_spinor_xmy(z, bprime, y);
    }
    else
    {
      /* Recursive update */
      alpha.re = -CDEVR(gamma, omega);
      alpha.im = -CDEVI(gamma, omega);
      qpb_spinor_axpy(r, alpha, w, r);   /* r -= alpha · D_toy · p  */
      qpb_spinor_axpy(z, alpha, y, z);   /* z -= alpha · D†D · p    */
    }

    /* res_norm = ||z_{k+1}||^2  (real, positive) */
    qpb_spinor_xdotx(&res_norm, z);

    /* beta = ||z_{k+1}||^2 / ||z_k||^2  (real) */
    beta.re = res_norm / gamma.re;
    beta.im = 0.;

    /* p_{k+1} = z_{k+1} + beta · p_k */
    qpb_spinor_axpy(p, beta, p, z);

    gamma.re = res_norm;
    gamma.im = 0.;

    qpb_spinor_xdotx(&true_res_norm, r);
  }

  if(iters == max_iter)
    return -1;

  return iters;
}


/* ================================================================
   Outer solver: preconditioned CGNE for  D_toy · x = b
   Preconditioner:  same D_toy, solved by qpb_toy_inner_CG.

   In this toy case the preconditioner is identical to the outer
   operator, so M^{-1}A = I.  The outer loop should converge in
   O(1) iterations with a sufficiently tight inner tolerance.

   Prints progress every n_echo iterations and final residual.
   Returns iteration count or -1 on non-convergence.

   Temp slots used: [7..13]
     [7]  r      residual of original system  b - D_toy·x
     [8]  p      search direction
     [9]  z      normal-equations residual    D†_toy · r
     [10] y      D†_toy · w  = (D†D) · p
     [11] w      D_toy · p
     [12] bprime D†_toy · b
     [13] s      output of inner preconditioner solve: M·s = z
   ================================================================ */
int
qpb_congrad_overlap_kl_pfrac(qpb_spinor_field x, qpb_spinor_field b,
                             qpb_double CG_epsilon,   int CG_max_iter)
{
  qpb_spinor_field r      = toy_temp_vecs[7];
  qpb_spinor_field p      = toy_temp_vecs[8];
  qpb_spinor_field z      = toy_temp_vecs[9];
  qpb_spinor_field y      = toy_temp_vecs[10];
  qpb_spinor_field w      = toy_temp_vecs[11];
  qpb_spinor_field bprime = toy_temp_vecs[12];
  qpb_spinor_field s      = toy_temp_vecs[13];

  int n_reeval = 100;
  int n_echo   = 10;     /* print every 10 iters so slow convergence is visible */
  int iters    = 0;

  qpb_double res_norm, true_res_norm, b_norm, bprime_norm;
  qpb_complex_double alpha = {1, 0}, omega = {1, 0};
  qpb_complex_double beta, gamma, new_gamma;

  /* ||b|| */
  qpb_spinor_xdotx(&b_norm, b);
  true_res_norm = b_norm;

  /* bprime = D†_toy · b */
  Dconj_toy_op(bprime, b);
  qpb_spinor_xdotx(&bprime_norm, bprime);

  /* x0 = 0 */
  qpb_spinor_field_set_zero(x);

  /* r0 = b,  z0 = bprime  (exact since x0 = 0) */
  qpb_spinor_xeqy(r, b);
  qpb_spinor_xeqy(z, bprime);

  /* Solve M·s0 = z0 for s0  (first preconditioner application) */
  int prec_iters = qpb_toy_inner_CG(s, z, prec_epsilon, prec_max_iter);
  if(prec_iters < 0)
    error(" [toy] WARNING: inner CG did not converge at initialization\n");

  /* gamma_0 = z0† · s0  (real by HPD of M = D†D; enforce .im = 0) */
  qpb_spinor_xdoty(&gamma, z, s);
  gamma.im = 0.;

  /* p0 = s0  (preconditioned initial search direction) */
  qpb_spinor_xeqy(p, s);

  print(" [toy] Starting preconditioned CGNE\n");
  print(" [toy]   outer tol = %e,  inner tol = %e\n", CG_epsilon, prec_epsilon);
  print(" [toy]   gamma_0 = %e  (should be positive and real)\n", gamma.re);

  qpb_double t = qpb_stop_watch(0);
  for(iters = 1; iters < CG_max_iter; iters++)
  {
    /* Stopping criterion on true residual of original system D_toy·x = b */
    if(true_res_norm / b_norm <= CG_epsilon)
      break;

    /* w = D_toy · p,  y = D†_toy · w = (D†D) · p */
    D_toy_op(w, p);
    Dconj_toy_op(y, w);

    /* omega = p† · (D†D) · p = ||D·p||^2 = ||w||^2  (real, positive) */
    // qpb_spinor_xdotx(&omega.re, w);
    // omega.im = 0.;
    qpb_spinor_xdoty(&omega, p, y);

    /* alpha = gamma / omega */
    alpha = CDEV(gamma, omega);

    /* x += alpha · p */
    qpb_spinor_axpy(x, alpha, p, x);

    /* Update r and z */
    if(iters % n_reeval == 0)
    {
      /* Full recomputation to suppress round-off */
      D_toy_op(w, x);
      qpb_spinor_xmy(r, b, w);
      Dconj_toy_op(y, w);
      qpb_spinor_xmy(z, bprime, y);
    }
    else
    {
      alpha.re = -CDEVR(gamma, omega);
      alpha.im = -CDEVI(gamma, omega);
      qpb_spinor_axpy(r, alpha, w, r);   /* r -= alpha · D_toy · p    */
      qpb_spinor_axpy(z, alpha, y, z);   /* z -= alpha · D†D_toy · p  */
    }

    /* Solve M·s = z  (one inner solve per outer iteration) */
    prec_iters = qpb_toy_inner_CG(s, z, prec_epsilon, prec_max_iter);
    if(prec_iters < 0)
      error(" [toy] WARNING: inner CG did not converge at outer iter %d\n", iters);

    /* new_gamma = z† · s  (real by HPD of M; take .re explicitly) */
    qpb_spinor_xdoty(&new_gamma, z, s);
    res_norm = new_gamma.re;

    /* beta = gamma_{k+1} / gamma_k  (real) */
    beta.re = res_norm / gamma.re;
    beta.im = 0.;

    /* p_{k+1} = s_{k+1} + beta · p_k */
    qpb_spinor_axpy(p, beta, p, s);

    /* Advance gamma */
    gamma.re = res_norm;
    gamma.im = 0.;

    qpb_spinor_xdotx(&true_res_norm, r);
    if(iters % n_echo == 0)
      print(" [toy] \t iters = %4d, res = %e, inner_iters = %d\n",
             iters, true_res_norm / b_norm, prec_iters);
  }
  t = qpb_stop_watch(t);

  /* Final explicit residual of original system */
  D_toy_op(w, x);
  qpb_spinor_xmy(r, b, w);
  qpb_spinor_xdotx(&true_res_norm, r);

  if(iters == CG_max_iter)
  {
    error(" [toy] !\n");
    error(" [toy] Preconditioned CGNE *did not* converge after %d iterations\n", iters);
    error(" [toy] residual = %e,  relative = %e,  t = %g sec\n",
           true_res_norm, true_res_norm / b_norm, t);
    error(" [toy] !\n");
    return -1;
  }

  print(" [toy] Converged after %d outer iters, res = %e, relative = %e, t = %g sec\n",
         iters, true_res_norm, true_res_norm / b_norm, t);

  return iters;
}