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
#include <stdlib.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

/*
  ov_temp_vecs layout (25 vectors total):

   Apply scratch (transient; clobbered on every call to the apply functions):
    [ 0]  sum vector for the partial-fraction accumulation (also reused as the
          low-mode scratch inside apply_gamma5_sign, after D_op has consumed it)
    [ 1]  z vector  (D·sum in apply_overlap)
    [ 2]  γ5 scratch in apply_conj_overlap

   Outer BiCGStab persistent state:
    [ 3]  r0_hat                         [ 7]  v = D_ov^outer · z_pc
    [ 4]  r                              [ 8]  y_pc  = K^{-1} · p
    [ 5]  p                              [ 9]  z_pc  = K^{-1} · s
    [ 6]  u = D_ov^outer · y_pc

   Preconditioner-solver scratch (one set of slots):
     n_outer → inner BiCGStab on D_ov^(n_outer-1):
       [10]  r0_hat_in   [11] r_in   [12] p_in   [13] u_in   [14] v_in
       [15] prec1: y_pc   ← NEW (K2⁻¹·p) [16] prec1: z_pc   ← NEW (K2⁻¹·s)

    [17] prec2: r0_hat   [18] prec2: r   [19] prec2: p   [20] prec2: u   [21] prec2: v

   Sign-function deflation scratch (transient apply scratch, like [0]/[1]/[2]):
    [24]  P_high·x = x - V V^dag x   (the deflated MSCG right-hand side)

*/

#define OVERLAP_NUMB_TEMP_VECS 25
#define MSCG_NUMB_TEMP_VECS    20

/* ---- Experimental level-2 cascaded preconditioning (hard-coded knobs) -------
   Deliberately compile-time constants, NOT parsed from the input file: flip a
   value and recompile. Experimental phase — may be removed entirely.            */
#define SECOND_LAYER_REQUESTED   0      /* 1 = enable L2, 0 = disable           */
#define PREC2_MS_EPSILON_FACTOR  5.0    /* L2 MSCG tol = factor x L1 MSCG tol   */
#define PREC2_MAX_ITER_OFFSET    1      /* L2 BiCGStab cap = L1 cap - offset    */
#define PREC2_EPSILON_FACTOR     5.0    /* L2 BiCGStab tol = factor x L1 BiCGStab tol */

/* ========================================================================= *
 *  Sign-function low-mode deflation  (mirrors lib/qpb_overlap_Zolotarev.c)
 * ========================================================================= *
 *
 *  The overlap sign function sign(X), X = gamma5 (D - rho), is hard because of
 *  the smallest eigenvalues of X (equivalently, of X^2): the partial-fraction
 *  rational R(X^2) ~ (X^2)^{-1/2} must resolve them, which is what makes the
 *  inner multi-shift CG expensive.  Exactly as in the Zolotarev module, we
 *  deflate those low modes OUT OF THE OPERATOR (not the linear solve):
 *
 *      sign(X) x = V sign(W) (V^dag x)          [low modes, treated exactly]
 *                + X R(X^2) (x - V V^dag x)      [rational on the complement]
 *
 *  where V holds the k lowest eigenvectors of X^2 (an X-invariant subspace),
 *  and W = V^dag X V is the k x k Hermitian projection of the (signed) kernel
 *  onto that subspace.  sign(W) = Q sign(Theta) Q^dag is computed from the dense
 *  eigendecomposition of W (robust to the +/-lambda degeneracy of X^2).
 *
 *  Only the multi-shift RIGHT-HAND SIDE is deflated (P_high x has no low-mode
 *  content); the multi-shift CG and its zeta-recurrence are untouched, so every
 *  partial-fraction pole still converges to its full tolerance.  The subspace
 *  depends only on X, so ONE build serves every KL level (outer, prec, prec2).
 *
 *  This is verbatim the Zolotarev construction; the ONLY differences here are:
 *    (i)  the KL init has no pre-existing Lanczos pass, so the build runs as a
 *         standalone pass (it needs no lambda_min/lambda_max), and
 *    (ii) the order-0 fallback (sign(X) ~ X, a bare D_op) has no MSCG and is
 *         left untouched — deflation applies only when kl_order[lvl] > 0.
 *
 *  Compile-time controls (override e.g. with -DQPB_DEFL_K=0 for the baseline):
 *      QPB_DEFL_K : number of low modes carried in V (0 disables deflation).
 *      QPB_DEFL_M : Lanczos search dimension used to build V (>= K).
 * ========================================================================= */
#ifndef QPB_DEFL_K
#define QPB_DEFL_K 20
#endif
#ifndef QPB_DEFL_M
#define QPB_DEFL_M 64
#endif
#define QPB_DEFL_KMAX (QPB_DEFL_K > 0 ? QPB_DEFL_K : 1)

typedef qpb_complex_double cdbl;

typedef enum {
  LEVEL_OUTER = 0,
  LEVEL_PREC  = 1,
  LEVEL_PREC2 = 2,
  LEVEL_COUNT = 3
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

/* Preconditioner-solver control: drives inner BiCGStab on D_ov^(n_outer-1). */
static qpb_double prec_solver_epsilon;
static int        prec_solver_max_iter;
static int        prec2_solver_max_iter;
static qpb_double prec2_solver_epsilon;

static int        second_layer_on;        /* resolved once at init */

/* Sign-function deflation subspace (built once per configuration). */
static qpb_spinor_field defl_V[QPB_DEFL_KMAX];                  /* low eigvecs of X^2 */
static cdbl             defl_signW[QPB_DEFL_KMAX*QPB_DEFL_KMAX]; /* sign(V^dag X V)   */
static int              defl_k = 0;          /* usable modes (0 => disabled) */
static int              defl_built = 0;      /* has V been built for this config? */


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


INLINE void
X_op(qpb_spinor_field y, qpb_spinor_field x)
{
  /* Hermitian Wilson/Brillouin kernel  X = gamma5 (D - rho).  Same operator
     the inner MSCG squares; used here to build the deflation subspace and its
     projected sign. */

  void *dslash_args[4];

  dslash_args[0] = ov_params.gauge_ptr;
  dslash_args[1] = &ov_params.m_bare;
  dslash_args[2] = &ov_params.clover;
  dslash_args[3] = &ov_params.c_sw;

  ov_params.g5_dslash_op(y, x, dslash_args);

  return;
}


/* ------------------------- DEFLATION SUBSPACE BUILD ---------------------- *
 *  Build V (the k lowest eigenvectors of A = X^2) by a Lanczos pass with full
 *  re-orthogonalization, lifting the k lowest Ritz vectors and orthonormalizing
 *  them; then form W = V^dag X V (k x k Hermitian) and its matrix sign,
 *  sign(W) = Q sign(Theta) Q^dag.  All persistent state (defl_V, defl_signW)
 *  is filled here, once per gauge configuration.  Identical to the Zolotarev
 *  module's build (the KL init has no separate extreme-eigenvalue Lanczos to
 *  fold into, so this stands alone; it needs no lambda_min/lambda_max).
 * ------------------------------------------------------------------------- */
static void
build_deflation_subspace(void)
{
  int k = QPB_DEFL_K;
  int m = QPB_DEFL_M;
  if(k <= 0)
    return;
  if(m < k)
    m = k;

  print(" Deflation: building %d low modes of X^2 via Lanczos (m=%d)...\n", k, m);
  qpb_double tb = qpb_stop_watch(0);

  /* Lanczos vectors + scratch */
  qpb_spinor_field *lv = qpb_alloc(sizeof(qpb_spinor_field)*m);
  for(int i=0; i<m; i++)
    {
      lv[i] = qpb_spinor_field_init();
      qpb_spinor_field_set_zero(lv[i]);
    }
  qpb_spinor_field av  = qpb_spinor_field_init();
  qpb_spinor_field tmp = qpb_spinor_field_init();
  qpb_spinor_field_set_zero(av);
  qpb_spinor_field_set_zero(tmp);

  double *alpha = qpb_alloc(sizeof(double)*m);
  double *beta  = qpb_alloc(sizeof(double)*m);

  /* v_0 = normalized random vector */
  qpb_spinor_field_set_random(lv[0]);
  qpb_double nrm;
  qpb_spinor_xdotx(&nrm, lv[0]);
  qpb_spinor_ax(lv[0], (qpb_complex){1./sqrt(nrm), 0.}, lv[0]);

  for(int i=0; i<m; i++)
    {
      /* av = A v_i = X (X v_i) */
      X_op(tmp, lv[i]);
      X_op(av,  tmp);

      if(i > 0)
	qpb_spinor_axpy(av, (qpb_complex){-beta[i-1], 0.}, lv[i-1], av);

      qpb_complex_double a;
      qpb_spinor_xdoty(&a, lv[i], av);
      alpha[i] = a.re;
      qpb_spinor_axpy(av, (qpb_complex){-alpha[i], 0.}, lv[i], av);

      /* full re-orthogonalization (twice, for numerical stability) */
      for(int pass=0; pass<2; pass++)
	for(int j=0; j<=i; j++)
	  {
	    qpb_complex_double c;
	    qpb_spinor_xdoty(&c, lv[j], av);
	    qpb_spinor_axpy(av, (qpb_complex){-c.re, -c.im}, lv[j], av);
	  }

      qpb_double bb;
      qpb_spinor_xdotx(&bb, av);
      beta[i] = sqrt(bb);

      if(i < m-1)
	{
	  if(beta[i] < 1e-12)
	    {
	      /* invariant subspace found -- restart with a fresh random vector
		 (it gets re-orthogonalized at the next step) */
	      qpb_spinor_field_set_random(lv[i+1]);
	    }
	  else
	    {
	      qpb_spinor_ax(lv[i+1], (qpb_complex){1./beta[i], 0.}, av);
	    }
	}
    }

  /* diagonalize the symmetric tridiagonal T (alpha[0..m-1], beta[0..m-2]) */
  gsl_matrix *T = gsl_matrix_calloc(m, m);
  for(int i=0; i<m; i++)
    gsl_matrix_set(T, i, i, alpha[i]);
  for(int i=0; i<m-1; i++)
    {
      gsl_matrix_set(T, i, i+1, beta[i]);
      gsl_matrix_set(T, i+1, i, beta[i]);
    }
  gsl_vector *eval = gsl_vector_alloc(m);
  gsl_matrix *evec = gsl_matrix_alloc(m, m);
  gsl_eigen_symmv_workspace *ws = gsl_eigen_symmv_alloc(m);
  gsl_eigen_symmv(T, eval, evec, ws);
  gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_VAL_ASC);
  gsl_eigen_symmv_free(ws);

  /* lift the k lowest Ritz vectors:  V_c = sum_i evec[i][c] v_i */
  for(int c=0; c<k; c++)
    {
      defl_V[c] = qpb_spinor_field_init();
      qpb_spinor_field_set_zero(defl_V[c]);
      for(int i=0; i<m; i++)
	{
	  double sic = gsl_matrix_get(evec, i, c);
	  qpb_spinor_axpy(defl_V[c], (qpb_complex){sic, 0.}, lv[i], defl_V[c]);
	}
    }

  gsl_matrix_free(T);
  gsl_vector_free(eval);
  gsl_matrix_free(evec);

  /* modified Gram-Schmidt: enforce V^dag V = I */
  for(int a=0; a<k; a++)
    {
      for(int b=0; b<a; b++)
	{
	  qpb_complex_double c;
	  qpb_spinor_xdoty(&c, defl_V[b], defl_V[a]);
	  qpb_spinor_axpy(defl_V[a], (qpb_complex){-c.re, -c.im},
			  defl_V[b], defl_V[a]);
	}
      qpb_double n2;
      qpb_spinor_xdotx(&n2, defl_V[a]);
      qpb_spinor_ax(defl_V[a], (qpb_complex){1./sqrt(n2), 0.}, defl_V[a]);
    }

  /* Ritz residuals ||X^2 V_c - theta_c V_c|| : accuracy check for the subspace.
     If the higher modes' residuals are loose, increase QPB_DEFL_M. */
  print(" Deflation: X^2 low Ritz values / residuals (post-orthonormalization):\n");
  for(int c=0; c<k; c++)
    {
      X_op(tmp, defl_V[c]);
      X_op(av,  tmp);                        /* av = X^2 V_c */
      qpb_complex_double th;
      qpb_spinor_xdoty(&th, defl_V[c], av);  /* theta_c = <V_c, X^2 V_c> */
      qpb_spinor_axpy(av, (qpb_complex){-th.re, 0.}, defl_V[c], av);
      qpb_double rr;
      qpb_spinor_xdotx(&rr, av);
      print("   mode %2d: theta = %+e, ||res|| = %e\n", c, th.re, sqrt(rr));
    }

  /* W = V^dag X V  (k x k Hermitian).  One kernel apply per column. */
  gsl_matrix_complex *Wm = gsl_matrix_complex_alloc(k, k);
  for(int j=0; j<k; j++)
    {
      X_op(tmp, defl_V[j]);                  /* tmp = X V_j */
      for(int i=0; i<k; i++)
	{
	  qpb_complex_double wij;
	  qpb_spinor_xdoty(&wij, defl_V[i], tmp);
	  gsl_matrix_complex_set(Wm, i, j, gsl_complex_rect(wij.re, wij.im));
	}
    }

  /* W = Q Theta Q^dag  (Theta real, signed kernel eigenvalues in the subspace) */
  gsl_vector *th = gsl_vector_alloc(k);
  gsl_matrix_complex *Q = gsl_matrix_complex_alloc(k, k);
  gsl_eigen_hermv_workspace *wsh = gsl_eigen_hermv_alloc(k);
  gsl_eigen_hermv(Wm, th, Q, wsh);
  gsl_eigen_hermv_free(wsh);

  /* sign(W)[i][l] = sum_m Q[i][m] sign(theta_m) conj(Q[l][m]) */
  for(int i=0; i<k; i++)
    for(int l=0; l<k; l++)
      {
	double sr = 0., si = 0.;
	for(int mm=0; mm<k; mm++)
	  {
	    double sg = (gsl_vector_get(th, mm) >= 0.) ? 1. : -1.;
	    gsl_complex qim = gsl_matrix_complex_get(Q, i, mm);
	    gsl_complex qlm = gsl_matrix_complex_get(Q, l, mm);
	    /* qim * conj(qlm) */
	    double pr = GSL_REAL(qim)*GSL_REAL(qlm) + GSL_IMAG(qim)*GSL_IMAG(qlm);
	    double pi = GSL_IMAG(qim)*GSL_REAL(qlm) - GSL_REAL(qim)*GSL_IMAG(qlm);
	    sr += sg*pr;
	    si += sg*pi;
	  }
	defl_signW[i*k+l] = (cdbl){sr, si};
      }

  print(" Deflation: signed kernel eigenvalues (H_W) in the low subspace:\n");
  for(int mm=0; mm<k; mm++)
    print("   lambda[%2d] = %+e\n", mm, gsl_vector_get(th, mm));

  /* sanity: ||sign(W)^2 - I||_F  (should be ~0; sign is an involution) */
  {
    double dev = 0.;
    for(int i=0; i<k; i++)
      for(int l=0; l<k; l++)
	{
	  double sr = 0., si = 0.;
	  for(int mm=0; mm<k; mm++)
	    {
	      cdbl a = defl_signW[i*k+mm];
	      cdbl b = defl_signW[mm*k+l];
	      sr += a.re*b.re - a.im*b.im;
	      si += a.re*b.im + a.im*b.re;
	    }
	  double dr = sr - (i==l ? 1. : 0.);
	  dev += dr*dr + si*si;
	}
    print(" Deflation: ||sign(W)^2 - I||_F = %e\n", sqrt(dev));
  }

  gsl_matrix_complex_free(Wm);
  gsl_vector_free(th);
  gsl_matrix_complex_free(Q);

  for(int i=0; i<m; i++)
    qpb_spinor_field_finalize(lv[i]);
  free(lv);
  qpb_spinor_field_finalize(av);
  qpb_spinor_field_finalize(tmp);
  free(alpha);
  free(beta);

  defl_k = k;
  defl_built = 1;

  tb = qpb_stop_watch(tb);
  print(" Deflation: subspace ready (k=%d), build t = %g secs\n", k, tb);

  return;
}


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

  /* Levels: outer = kl_iters, prec = kl_iters - 1.
    At n_outer = 1, LEVEL_PREC has order 0 and apply_gamma5_sign uses
    sign(X) ≈ X  ⟹  D_ov^(0) = C·I + R·(D − ρ)  (shifted kernel). */
  kl_order[LEVEL_OUTER] = kl_iters;
  kl_order[LEVEL_PREC]  = kl_iters - 1;
  kl_order[LEVEL_PREC2] = kl_iters - 2;          /* may be < 0 for kl_iters < 2 */

  MS_solver_precision[LEVEL_OUTER] = ms_epsilon;
  MS_solver_precision[LEVEL_PREC]  = prec_ms_epsilon;
  MS_solver_precision[LEVEL_PREC2] = PREC2_MS_EPSILON_FACTOR * prec_ms_epsilon;

  MS_maximum_solver_iterations     = ms_max_iters;

  prec_solver_epsilon  = prec_epsilon;
  prec_solver_max_iter = prec_max_iter;

  prec2_solver_max_iter = prec_solver_max_iter - PREC2_MAX_ITER_OFFSET;
  prec2_solver_epsilon  = PREC2_EPSILON_FACTOR * prec_solver_epsilon;

  /* Populate partial-fraction tables per level. Skip LEVEL_PREC when its
     order is 0 — the kernel path doesn't use shifts/numerators. */
  for(int lvl = 0; lvl < LEVEL_COUNT; lvl++) {
    int n = kl_order[lvl];
    if(n <= 0) { shifts[lvl] = NULL; numerators[lvl] = NULL;
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

  /* Second layer is meaningful only if the first layer is actually running. */
  second_layer_on = SECOND_LAYER_REQUESTED
                    && (prec_max_iter > 0)
                    && (kl_order[LEVEL_PREC2] >= 0)
                    && (prec2_solver_max_iter > 1);

  /* ------------------- sign-function deflation subspace ------------------ *
     Build V (lowest k eigenvectors of X^2) and sign(V^dag X V) once, here.
     Reused by apply_gamma5_sign() at every level (they share X).  Disabled at
     compile time with -DQPB_DEFL_K=0, which reproduces the plain path. */
  if(QPB_DEFL_K > 0)
    build_deflation_subspace();

  /* MSCG workspace sized for the larger of the KL orders. */
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

  /* release the deflation subspace */
  if(defl_built) {
    for(int i=0; i<defl_k; i++)
      qpb_spinor_field_finalize(defl_V[i]);
    defl_built = 0;
    defl_k = 0;
  }
}


static void
apply_gamma5_sign(overlap_level_t lvl,
                  qpb_spinor_field y, qpb_spinor_field x)
{
  /* γ5·sign(X)·x at the partial-fraction order of this level.

     With low-mode deflation (defl_built, k > 0):

        γ5 sign(X) x =  D_op( R(X^2) P_high x )         [ high modes, rational ]
                     +  γ5 V sign(W) (V^dag x)          [ low modes, exact     ]

     where P_high x = x - V (V^dag x).  Only the RHS is deflated, so the MSCG
     zeta-recurrence is untouched and every pole keeps its full accuracy. */

  int n = kl_order[lvl];
  if(n == 0) {
    /* sign(X) ≈ X  ⇒  γ5·sign(X)·x = γ5·X·x = (D − ρ)·x .
       No MSCG here, hence nothing to deflate — left exactly as before. */
    D_op(y, x);
    return;
  }

  int k = defl_built ? defl_k : 0;

  qpb_spinor_field sum = ov_temp_vecs[0];
  /* MSCG right-hand side: the deflated P_high x when k>0, else x itself. */
  qpb_spinor_field rhs = (k > 0) ? ov_temp_vecs[24] : x;

  qpb_spinor_field yMS[n];
  for(int s = 0; s < n; s++) {
    yMS[s] = mscg_temp_vecs[s];
    qpb_spinor_field_set_zero(yMS[s]);
  }

  /* ---- (a) deflate the right-hand side:  a = V^dag x ; rhs = x - V a ---- */
  cdbl a[QPB_DEFL_KMAX];
  if(k > 0) {
    for(int i = 0; i < k; i++)
      qpb_spinor_xdoty(&a[i], defl_V[i], x);        /* a_i = <V_i, x> */
    qpb_spinor_xeqy(rhs, x);
    for(int i = 0; i < k; i++)
      qpb_spinor_axpy(rhs, (qpb_complex){-a[i].re, -a[i].im}, defl_V[i], rhs);
  }

  qpb_double kernel_mass  = ov_params.m_bare;
  qpb_double kernel_kappa = 1.0 / (2*kernel_mass + 8.0);

  qpb_mscongrad(yMS, rhs, ov_params.gauge_ptr, ov_params.clover, kernel_kappa,
                n, shifts[lvl], ov_params.c_sw,
                MS_solver_precision[lvl], MS_maximum_solver_iterations);

  /* ---- (b) rational part on the (deflated) RHS:  sum = R(X^2) rhs ---- */
  qpb_spinor_ax(sum, (qpb_complex){constant_term[lvl], 0.}, rhs);
  for(int s = 0; s < n; s++)
    qpb_spinor_axpy(sum, (qpb_complex){numerators[lvl][s], 0.},
                    yMS[s], sum);

  D_op(y, sum);   /* high-mode result: γ5 X R(X^2) P_high x */

  /* ---- (c) exact low-mode part:  y += γ5 V sign(W) a ---- */
  if(k > 0) {
    cdbl b[QPB_DEFL_KMAX];
    for(int i = 0; i < k; i++) {
      cdbl acc = {0., 0.};
      for(int j = 0; j < k; j++) {
        cdbl w = defl_signW[i*k+j];
        acc.re += w.re*a[j].re - w.im*a[j].im;
        acc.im += w.re*a[j].im + w.im*a[j].re;
      }
      b[i] = acc;
    }
    /* sum ([0]) is dead after D_op -- reuse it to accumulate the low part */
    qpb_spinor_field low = ov_temp_vecs[0];
    qpb_spinor_field_set_zero(low);
    for(int i = 0; i < k; i++)
      qpb_spinor_axpy(low, (qpb_complex){b[i].re, b[i].im}, defl_V[i], low);
    qpb_spinor_gamma5(low, low);
    qpb_spinor_axpy(y, (qpb_complex){1., 0.}, low, y);   /* y += γ5 V sign(W) a */
  }
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
preconditioner_bicgstab_2(qpb_spinor_field x, qpb_spinor_field b)
{
  qpb_spinor_field r0_hat = ov_temp_vecs[17];
  qpb_spinor_field r      = ov_temp_vecs[18];
  qpb_spinor_field p      = ov_temp_vecs[19];
  qpb_spinor_field u      = ov_temp_vecs[20];
  qpb_spinor_field v      = ov_temp_vecs[21];

  qpb_double res_norm, b_norm;
  qpb_complex_double alpha = {1, 0}, omega = {1, 0};
  qpb_complex_double beta, gamma, rho, zeta;
  int iters;

  qpb_spinor_xdotx(&b_norm, b);

  qpb_spinor_field_set_zero(x);
  qpb_spinor_field_set_zero(p);
  qpb_spinor_field_set_zero(u);

  qpb_spinor_xeqy(r, b);            /* r0 = b - D_ov^prec·x0 = b */
  qpb_spinor_xeqy(r0_hat, r);

  qpb_spinor_xdotx(&res_norm, r);
  gamma.re = res_norm; gamma.im = 0.;
  rho = gamma;

  for(iters = 1; iters < prec2_solver_max_iter; iters++) {
    if(res_norm / b_norm <= prec2_solver_epsilon)
      break;

    qpb_spinor_xdoty(&gamma, r0_hat, r);
    beta = CMUL(CDEV(gamma, rho), CDEV(alpha, omega));

    omega = CNEGATE(omega);
    qpb_spinor_axpy(p, omega, u, p);     /* p -= omega·u                 */
    qpb_spinor_axpy(p, beta, p, r);      /* p  = beta·p + r              */

    apply_overlap(LEVEL_PREC2, u, p);     /* u  = D_ov^prec · p           */

    qpb_spinor_xdoty(&beta, r0_hat, u);
    rho   = gamma;
    alpha = CDEV(rho, beta);

    alpha = CNEGATE(alpha);
    qpb_spinor_axpy(r, alpha, u, r);     /* r -= alpha·u   (r ≡ s now)   */

    apply_overlap(LEVEL_PREC2, v, r);     /* v  = D_ov^prec · s           */

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

  print(" \t\t\tpreconditioner-2 BiCGStab: %d iters, relative residual = %e\n",
        iters, res_norm / b_norm);

  return iters;            /* caller can compare against prec_solver_max_iter */
}

static int
preconditioner_bicgstab(qpb_spinor_field x, qpb_spinor_field b)
{
  qpb_spinor_field r0_hat = ov_temp_vecs[10];
  qpb_spinor_field r      = ov_temp_vecs[11];
  qpb_spinor_field p      = ov_temp_vecs[12];
  qpb_spinor_field u      = ov_temp_vecs[13];
  qpb_spinor_field v      = ov_temp_vecs[14];
  qpb_spinor_field y_pc   = ov_temp_vecs[15];
  qpb_spinor_field z_pc   = ov_temp_vecs[16];

  qpb_double res_norm, b_norm;
  qpb_complex_double alpha = {1, 0}, omega = {1, 0};
  qpb_complex_double beta, gamma, rho, zeta;
  int iters;

  qpb_spinor_xdotx(&b_norm, b);

  qpb_spinor_field_set_zero(x);
  qpb_spinor_field_set_zero(p);
  qpb_spinor_field_set_zero(u);

  qpb_spinor_xeqy(r, b);            /* r0 = b - D_ov^prec·x0 = b */
  qpb_spinor_xeqy(r0_hat, r);

  qpb_spinor_xdotx(&res_norm, r);
  gamma.re = res_norm; gamma.im = 0.;
  rho = gamma;

  for(iters = 1; iters < prec_solver_max_iter; iters++) {
    if(res_norm / b_norm <= prec_solver_epsilon)
      break;

    qpb_spinor_xdoty(&gamma, r0_hat, r);
    beta = CMUL(CDEV(gamma, rho), CDEV(alpha, omega));

    omega = CNEGATE(omega);
    qpb_spinor_axpy(p, omega, u, p);     /* p -= omega·u                 */
    qpb_spinor_axpy(p, beta, p, r);      /* p  = beta·p + r              */

    // apply_overlap(LEVEL_PREC, u, p);     /* u  = D_ov^prec · p           */
    if(second_layer_on) preconditioner_bicgstab_2(y_pc, p);  /* y_pc = K2⁻¹·p   */
    else                qpb_spinor_xeqy(y_pc, p);
    apply_overlap(LEVEL_PREC, u, y_pc);                       /* u = D_ov^(n-1)·y_pc */

    qpb_spinor_xdoty(&beta, r0_hat, u);
    rho   = gamma;
    alpha = CDEV(rho, beta);

    alpha = CNEGATE(alpha);
    qpb_spinor_axpy(r, alpha, u, r);     /* r -= alpha·u   (r ≡ s now)   */

    // apply_overlap(LEVEL_PREC, v, r);     /* v  = D_ov^prec · s           */
    if(second_layer_on) preconditioner_bicgstab_2(z_pc, r);  /* z_pc = K2⁻¹·s   */
    else                qpb_spinor_xeqy(z_pc, r);
    apply_overlap(LEVEL_PREC, v, z_pc);                       /* v = D_ov^(n-1)·z_pc */

    qpb_spinor_xdoty(&zeta, v, r);
    qpb_spinor_xdotx(&beta.re, v); beta.im = 0;
    omega = CDEV(zeta, beta);

    alpha = CNEGATE(alpha);
    // qpb_spinor_axpy(x, alpha, p, x);     /* x += alpha·p                 */
    // qpb_spinor_axpy(x, omega, r, x);     /* x += omega·s   (s is in r)   */
    qpb_spinor_axpy(x, alpha, y_pc, x);    /* x += alpha·y_pc   (NOT alpha·p) */
    qpb_spinor_axpy(x, omega, z_pc, x);    /* x += omega·z_pc   (NOT omega·s) */

    omega = CNEGATE(omega);
    qpb_spinor_axpy(r, omega, v, r);     /* r -= omega·v                 */
    omega = CNEGATE(omega);

    qpb_spinor_xdotx(&res_norm, r);
  }

  /* Explicit final (true) residual, reported for diagnostics — the
     recurrence residual used in the exit test can drift from b - D_ov^prec·x. */
  apply_overlap(LEVEL_PREC, u, x);
  qpb_spinor_xmy(r, b, u);
  qpb_spinor_xdotx(&res_norm, r);
  print(" \t\tpreconditioner BiCGStab: %d iters, relative residual = %e\n",
        iters, res_norm / b_norm);

  return iters;            /* caller can compare against prec_solver_max_iter */
}



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
  enum { PREC_NONE, PREC_BICGSTAB } prec_path;
  prec_path = (prec_solver_max_iter == 0) ? PREC_NONE : PREC_BICGSTAB;

  qpb_spinor_xdotx(&b_norm, b);
  qpb_spinor_field_set_zero(x);
  qpb_spinor_field_set_zero(p);
  qpb_spinor_field_set_zero(u);

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
