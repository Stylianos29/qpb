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
#include <qpb.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_sf_ellint.h>
#include <gsl/gsl_sf_elljac.h>


/*
  ov_temp_vecs layout (27 vectors total):

   Apply scratch (transient; clobbered on every call to the apply functions,
   shared across all levels since apply_overlap is atomic):
    [ 0]  sum vector for the partial-fraction accumulation (also reused as the
          low-mode scratch inside apply_gamma5_sign, after D_op has consumed it)
    [ 1]  z vector (D·sum in apply_overlap)

   CG-on-normal-equations solver state (qpb_congrad_overlap_Zolotarev, unchanged):
    [ 2]  p   [ 3]  r   [ 4]  y   [ 5]  w   [ 6]  bprime

   Outer right-preconditioned BiCGStab state (qpb_bicgstab_overlap_Zolotarev):
    [ 7]  r0_hat   [ 8]  r   [ 9]  p   [10]  u   [11]  v
    [12]  y_pc = K^{-1}·p    [13]  z_pc = K^{-1}·s

   L1 preconditioner BiCGStab state on D_ov^(n-1):
    [14]  r0_hat   [15]  r   [16]  p   [17]  u   [18]  v
    [19]  y_pc = K2^{-1}·p   [20]  z_pc = K2^{-1}·s
          (L2 hooks; when the second layer is OFF, y_pc≡p and z_pc≡s, so L1 is
           behaviourally identical to the single-layer version)

   L2 preconditioner BiCGStab state on D_ov^(n-2):
    [21]  r0_hat   [22]  r   [23]  p   [24]  u   [25]  v

   Sign-function deflation scratch (transient apply scratch, like [0]/[1]):
    [26]  P_high·x = x - V V^dag x   (the deflated MSCG right-hand side)

   The CG block [2-6] and the BiCGStab blocks [7-25] are never live at the same
   time (mutually-exclusive top-level entry points), kept disjoint for clarity.
*/

#define OVERLAP_NUMB_TEMP_VECS 27
#define MSCG_NUMB_TEMP_VECS 20

/* ---- Experimental level-2 cascaded preconditioning (hard-coded knobs) -------
   Deliberately compile-time constants, NOT parsed from the input file: flip a
   value and recompile. Mirrors lib/qpb_overlap_kl_pfrac.c. Experimental.        */
#define SECOND_LAYER_REQUESTED   0      /* 1 = enable L2, 0 = disable           */
#define PREC2_MS_EPSILON_FACTOR  5.0    /* L2 MSCG tol = factor x L1 MSCG tol   */
#define PREC2_MAX_ITER_OFFSET    1      /* L2 BiCGStab cap = L1 cap - offset    */


/* ========================================================================= *
 *  Sign-function low-mode deflation
 * ========================================================================= *
 *
 *  The overlap sign function sign(X), X = gamma5 (D - rho), is hard precisely
 *  because of the smallest eigenvalues of X (equivalently, of X^2): the
 *  Zolotarev rational R(X^2) ~ (X^2)^{-1/2} must resolve them, which is what
 *  makes the inner multi-shift CG expensive and forces a high rational order.
 *
 *  We deflate those low modes OUT OF THE OPERATOR (not the linear solve):
 *
 *      sign(X) x = V sign(W) (V^dag x)          [low modes, treated exactly]
 *                + X R(X^2) (x - V V^dag x)      [rational on the complement]
 *
 *  where V holds the k lowest eigenvectors of X^2 (an X-invariant subspace),
 *  and W = V^dag X V is the k x k Hermitian projection of the (signed) kernel
 *  onto that subspace.  sign(W) = Q sign(Theta) Q^dag is computed from the
 *  dense eigendecomposition of W (robust to the +/-lambda degeneracy of X^2).
 *
 *  Because we clean only the RIGHT-HAND SIDE (P_high x has no low-mode
 *  content), the multi-shift CG and its zeta-recurrence are untouched: every
 *  Zolotarev pole still converges to its full tolerance (no residual floor).
 *  The k eigenvectors are built once per gauge configuration and reused across
 *  every sign evaluation of ALL levels (outer, L1 prec, L2 prec), since they
 *  depend only on X.
 *
 *  When deflation is on (QPB_DEFL_K > 0) the build's Lanczos tridiagonal also
 *  supplies the extreme eigenvalues of X^2 for the Zolotarev interval, so the
 *  separate qpb_extreme_eigenvalues_of_X_squared() pass is skipped entirely.
 *
 *  Compile-time controls (override e.g. with -DQPB_DEFL_K=0 for the baseline):
 *      QPB_DEFL_K : number of low modes carried in V (0 disables deflation).
 *      QPB_DEFL_M : Lanczos search dimension used to build V (>= K).
 *      QPB_DEFL_SHRINK_INTERVAL : placeholder ONLY (default 0). Reserved for a
 *          future option to use [theta_k, lambda_max^2] as the Zolotarev
 *          interval once the low modes are handled exactly. It is intentionally
 *          NOT wired to anything yet: the interval is unchanged.
 * ========================================================================= */
#ifndef QPB_DEFL_K
#define QPB_DEFL_K 20
#endif
#ifndef QPB_DEFL_M
#define QPB_DEFL_M 64
#endif
#ifndef QPB_DEFL_SHRINK_INTERVAL
#define QPB_DEFL_SHRINK_INTERVAL 0
#endif
#define QPB_DEFL_KMAX (QPB_DEFL_K > 0 ? QPB_DEFL_K : 1)

typedef qpb_complex_double cdbl;


typedef enum {
  LEVEL_OUTER = 0,                 /* Zolotarev order n   (massive overlap)     */
  LEVEL_PREC  = 1,                 /* Zolotarev order n-1 (L1 preconditioner)   */
  LEVEL_PREC2 = 2,                 /* Zolotarev order n-2 (L2 preconditioner)   */
  LEVEL_COUNT = 3
} overlap_level_t;


static qpb_spinor_field ov_temp_vecs[OVERLAP_NUMB_TEMP_VECS];
static qpb_spinor_field mscg_temp_vecs[MSCG_NUMB_TEMP_VECS];

static qpb_overlap_params ov_params;

/* Per-level Zolotarev partial-fraction data. All levels share the SAME
   spectral interval [min_eigv, max_eigv]: the extreme eigenvalues of X^2 are a
   property of the operator, not of the Zolotarev order, so the (expensive)
   Lanczos pass runs once and every coefficient table is derived from it.
   For an order-0 level, shifts/numerators are NULL and constant_term holds the
   minimax linear coefficient c0 = 2/(min_eigv+max_eigv) (sign(X) ~ c0 X). */
static int          Zolotarev_order    [LEVEL_COUNT];
static qpb_double  *shifts             [LEVEL_COUNT];
static qpb_double  *numerators         [LEVEL_COUNT];
static qpb_double   constant_term      [LEVEL_COUNT];

/* MSCG control: per-level tolerance, shared iteration cap. */
static qpb_double   MS_solver_precision[LEVEL_COUNT];
static int          MS_maximum_solver_iterations;

/* Preconditioner-solver control. The inner BiCGStab iterates until its
   relative residual ||r||/||b|| <= prec_solver_epsilon, capped at
   prec_solver_max_iter steps (a hard safety net, NOT the stopping criterion).
   Because the iteration count depends on the right-hand side, the
   preconditioner is non-stationary; the outer plain BiCGStab tolerates this
   in practice provided prec_solver_epsilon is tight enough. */
static qpb_double   prec_solver_epsilon;     /* inner stopping criterion       */
static int          prec_solver_max_iter;    /* L1 hard safety cap             */
static int          prec2_solver_max_iter;   /* L2 hard safety cap             */
static int          prec_on;                 /* resolved once at init          */
static int          second_layer_on;         /* resolved once at init          */

/* Sign-function deflation subspace (built once per configuration). */
static qpb_spinor_field defl_V[QPB_DEFL_KMAX];                  /* low eigvecs of X^2 */
static cdbl             defl_signW[QPB_DEFL_KMAX*QPB_DEFL_KMAX]; /* sign(V^dag X V)   */
static int              defl_k = 0;          /* usable modes (0 => disabled) */
static int              defl_built = 0;      /* has V been built for this config? */


/* -------------- SCALAR FUNCTIONS -------------- */

/* Calculate c coefficients using Jacobi elliptic functions */
void
calculate_c_coefficients(qpb_double *c, int n, qpb_double k_prime, qpb_double K_prime)
{
    qpb_double sn, cn, dn;
    qpb_double u;
    qpb_double k2_prime = k_prime * k_prime;

    for (int i = 1; i <= 2*n; i++) {
        u = i * K_prime / (2.0 * n + 1.0);

        /* Compute Jacobi elliptic functions sn(u; k') */
        gsl_sf_elljac_e(u, k2_prime, &sn, &cn, &dn);

        qpb_double sn2 = sn * sn;
        c[i-1] = sn2 / (1.0 - sn2);
    }
}


/* Calculate b coefficients using c's */
void
calculate_b_coefficients(qpb_double *b, qpb_double *c, int n)
{
    for (int i = 1; i <= n; i++)
    {
      /* Numerator: product over k=1 to n of (c_{2k} - c_{2i-1}) */
      qpb_double numerator = 1.0;
      for (int k = 1; k <= n; k++)
        numerator *= c[2*k-1] - c[2*i-2];

      qpb_double denominator = 1.0;
      /* Denominator: product over k=1 to n (k≠i) of (c_{2k-1} - c_{2i-1}) */
      for (int k = 1; k <= n; k++)
        if (k != i)
          denominator *= c[2*k-2] - c[2*i-2];

      b[i-1] = numerator / denominator;
    }
}


qpb_double
sign_function_product_form(qpb_double x, int n, qpb_double *c)
{
  qpb_double x2 = x * x;
  qpb_double numerator = 1.0;
  qpb_double denominator = 1.0;

  for (int i = 0; i < n; i++)
  {
      numerator *= (x2 + c[2*i + 1]);
      denominator *= (x2 + c[2*i]);
  }

  return x * numerator / denominator;
}


// Calculate normalization constant
qpb_double
calculate_normalization_constant(qpb_double *c, int n, \
        qpb_double minimum_eigenvalue, qpb_double maximum_eigenvalue)
{
  qpb_double sign_function_at_1 = sign_function_product_form(1.0, n, c);
  qpb_double condition_number = maximum_eigenvalue / minimum_eigenvalue;
  qpb_double sign_function_at_max_range = sign_function_product_form(
                                                    condition_number, n, c);

  return 2.0 / (sign_function_at_1 + sign_function_at_max_range);
}


/* --------------------- EXTREME EIGENVALUES FUNCTIONS --------------------- */

INLINE void
tridiag_eigenv(qpb_double *eig, qpb_double *a, qpb_double *b, int n)
{
  /* It calculates the set of eigenvalues of the tri-diagonal matrix
  constructed appropriately from the given a and b arrays. */

  gsl_matrix *A = gsl_matrix_calloc(n, n);
  gsl_matrix_set (A, 0, 0, a[0]);
  gsl_matrix_set (A, 0, 0+1, b[0]);
  for(int i=1; i<n-1; i++)
    {
      gsl_matrix_set(A, i, i, a[i]);
      gsl_matrix_set(A, i, i+1, b[i]);
      gsl_matrix_set(A, i, i-1, b[i-1]);
    }
  gsl_matrix_set(A, n-1, n-1, a[n-1]);
  gsl_matrix_set(A, n-1, n-1-1, b[n-1-1]);

  gsl_vector *e = gsl_vector_alloc(n);
  gsl_eigen_symm_workspace *w = gsl_eigen_symm_alloc(n);
  gsl_eigen_symm(A, e, w);
  gsl_eigen_symm_free(w);
  gsl_matrix_free(A);

  gsl_sort_vector(e);

  for(int i=0; i<n; i++)
    eig[i] = gsl_vector_get(e, i);

  gsl_vector_free(e);

  return;
}


int
qpb_extreme_eigenvalues_of_X_squared(qpb_double *min_eigv, \
  qpb_double *max_eigv, qpb_double Lanczos_epsilon, int max_iters)
{
  /* It calculates the extreme eigenvalues of the eigenvalue spectrum
  of H^2, H ≡ γ5*Kernel(x), with: Kernel(x) = (a*D - ρ)(x), using the Lanczos
  algorithm.  Used only on the baseline (QPB_DEFL_K = 0) path; when deflation
  is on, the deflation build supplies these bounds instead. */

  qpb_lanczos_init();

  qpb_clover_term clover_term = ov_params.clover;
  qpb_double c_sw = ov_params.c_sw;
  qpb_double mass = ov_params.m_bare; // Kernel operator mass set at -rho
  qpb_double kappa = 1./(2*mass+8.);
  void *solver_arg_links = ov_params.gauge_ptr;

  qpb_double *a, *b, *eig;
  a = qpb_alloc(sizeof(qpb_double)*max_iters);
  b = qpb_alloc(sizeof(qpb_double)*max_iters);
  eig = qpb_alloc(sizeof(qpb_double)*max_iters);

  qpb_lanczos(a, b, solver_arg_links, clover_term, kappa, c_sw, 1);
  qpb_double lambda = 0, dlambda, lambda0 = 1e3;
  int i=0;
  for(i=1; i<max_iters; i++)
  {
    qpb_lanczos(a, b, solver_arg_links, clover_term, kappa, c_sw, -1);
    tridiag_eigenv(eig, a, b, i+1);

    lambda = eig[i] / eig[0];
    dlambda = fabs(lambda - lambda0) / fabs(lambda + lambda0);
    if (i%100==0)
      print("\titer = %4d, CN = %e/%e = %e (change = %e, target = %e)\n", i+1,\
                      eig[i], eig[0], eig[i]/eig[0], dlambda, Lanczos_epsilon);
    if(dlambda < Lanczos_epsilon*0.5)
      break;
    lambda0 = lambda;
  }

  *min_eigv = (qpb_double) eig[0];
  *max_eigv = (qpb_double) eig[i-1];

  return i;
}


/* ----------------------- EXPANSION COEFFICIENTS ----------------------- */

/* Build the Zolotarev partial-fraction tables (shifts, numerators,
   constant_term) for a given level/order, using the shared spectral interval
   [ov_params.min_eigv, ov_params.max_eigv]. Called once per active level with
   order >= 1. */
static void
build_zolotarev_tables(overlap_level_t lvl, int order, qpb_double scaling_factor)
{
  qpb_double min_eigv = ov_params.min_eigv;
  qpb_double max_eigv = ov_params.max_eigv;
  qpb_double min_eigv_squared = min_eigv * min_eigv;
  qpb_double max_eigv_squared = max_eigv * max_eigv;

  shifts    [lvl] = qpb_alloc(sizeof(qpb_double)*order);
  numerators[lvl] = qpb_alloc(sizeof(qpb_double)*order);

  /* Compute elliptic parameters from the shared spectral interval */
  qpb_double k_squared = min_eigv_squared/max_eigv_squared;
  qpb_double k_prime = sqrt(1.0 - k_squared);

  /* Compute complete elliptic integral */
  qpb_double K_prime = gsl_sf_ellint_Kcomp(k_prime, GSL_PREC_DOUBLE);

  /* Allocate scratch arrays */
  qpb_double *c = qpb_alloc(sizeof(qpb_double)*2*order);
  qpb_double *b = qpb_alloc(sizeof(qpb_double)*order);

  /* Compute coefficients */
  calculate_c_coefficients(c, order, k_prime, K_prime);
  calculate_b_coefficients(b, c, order);

  /* Compute normalization constant */
  qpb_double normalization_constant = calculate_normalization_constant(
                                c, order, min_eigv, max_eigv);

  constant_term[lvl] = normalization_constant / min_eigv;

  for(int i=0; i<order; i++)
  {
    shifts    [lvl][i] = c[2*i] * min_eigv_squared;
    numerators[lvl][i] = normalization_constant * b[i] * min_eigv;
  }

  /* Modify the numerical constants of the partial fraction expansion using the
     scaling parameter */
  if (scaling_factor != 1.0)
  {
    constant_term[lvl] *= 1/sqrt(scaling_factor);
    for(int i=0; i<order; i++)
    {
      numerators[lvl][i] *= sqrt(scaling_factor);
      shifts    [lvl][i] *= scaling_factor;
    }
  }

  free(c);
  free(b);
}


/* ------------------------ MATRIX-VECTOR FUNCTIONS ------------------------ */

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
 *  is filled here, once per gauge configuration.
 *
 *  Also returns the extreme eigenvalues of X^2 (edges of the Lanczos
 *  tridiagonal spectrum, which converge fastest), so the caller can use them
 *  as the Zolotarev spectral interval instead of a second, redundant pass.
 * ------------------------------------------------------------------------- */
static void
build_deflation_subspace(qpb_double *min_eigv_squared_out,
                         qpb_double *max_eigv_squared_out)
{
  int k = QPB_DEFL_K;
  int m = QPB_DEFL_M;
  if(k <= 0)
    return;
  if(m < k)
    m = k;

  print(" Deflation: building %d low modes of X^2 via Lanczos (m=%d)...\n", k, m);
  qpb_double tb = qpb_stop_watch(0);

  int kapps = 0;   /* count of forward kernel (X) applications in the build */

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
      kapps += 2;

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

  /* Harvest the spectral interval of X^2 from the tridiagonal edges. */
  if(min_eigv_squared_out) *min_eigv_squared_out = gsl_vector_get(eval, 0);
  if(max_eigv_squared_out) *max_eigv_squared_out = gsl_vector_get(eval, m-1);

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
      kapps += 2;
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
      kapps += 1;
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
  print(" Deflation: kernel applications for subspace build = %d\n", kapps);
  print(" Deflation: subspace ready (k=%d), build t = %g secs\n", k, tb);

  return;
}


void
qpb_overlap_Zolotarev_init(void * gauge, qpb_clover_term clover, \
          int Zol_iters, qpb_double rho, \
          qpb_double c_sw, qpb_double mass, qpb_double scaling_factor, \
          qpb_double ms_epsilon, qpb_double prec_ms_epsilon, int ms_max_iter, \
          qpb_double prec_epsilon, int prec_max_iter, \
          qpb_double Lanczos_epsilon, int Lanczos_max_iters, \
          qpb_double delta_max, qpb_double delta_min)
{
  if(ov_params.initialized != QPB_OVERLAP_INITIALIZED)
  {
    qpb_comm_halo_spinor_field_init();
    for(int i=0; i<OVERLAP_NUMB_TEMP_VECS; i++)
    {
      ov_temp_vecs[i] = qpb_spinor_field_init();
      qpb_spinor_field_set_zero(ov_temp_vecs[i]);
    }

    for(int i=0; i<MSCG_NUMB_TEMP_VECS; i++)
    {
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

    /* --------------------- spectral interval of X^2 --------------------- */

    qpb_double min_eigv_squared;
    qpb_double max_eigv_squared;

    /* The extrema of the eigenvalue spectrum of X^2, X = g5*(D - rho), feed
       the Zolotarev coefficient tables. They are a property of the operator,
       NOT of the Zolotarev order, so they are computed once.

       When deflation is on (QPB_DEFL_K > 0), build the subspace first and take
       the interval from its Lanczos tridiagonal edges (which converge fastest),
       avoiding a second, redundant extreme-eigenvalue pass. Otherwise use the
       standalone adaptive routine. */
    if(QPB_DEFL_K > 0)
    {
      build_deflation_subspace(&min_eigv_squared, &max_eigv_squared);
      print(" Total number of Lanczos algorithm iterations = %d\n", QPB_DEFL_M);
    }
    else
    {
      int Lanczos_iters = qpb_extreme_eigenvalues_of_X_squared(&min_eigv_squared,\
                        &max_eigv_squared, Lanczos_epsilon, Lanczos_max_iters);
      print(" Total number of Lanczos algorithm iterations = %d\n", \
                                                                Lanczos_iters);
    }

    /* If requested the extreme eigenvalues are modified accordingly */
    if (delta_min != 1.0)
      min_eigv_squared *= delta_min;
    if (delta_max != 1.0)
      max_eigv_squared *= delta_max;

    print(" Min eigenvalue squared = %.16f\n", min_eigv_squared);
    print(" Max eigenvalue squared = %.16f\n", max_eigv_squared);

    ov_params.min_eigv = sqrt(min_eigv_squared);
    ov_params.max_eigv = sqrt(max_eigv_squared);

    /* ----------------------- expansion coefficients ----------------------- */

    /* Levels: outer = Zol_iters, L1 = Zol_iters-1, L2 = Zol_iters-2.
       Order 0 is a valid member here: sign(X) ≈ c0 X (minimax linear), which
       lets an order-0 overlap act as the smallest preconditioner. */
    Zolotarev_order[LEVEL_OUTER] = Zol_iters;
    Zolotarev_order[LEVEL_PREC]  = Zol_iters - 1;
    Zolotarev_order[LEVEL_PREC2] = Zol_iters - 2;   /* may be < 0 for small n */

    MS_solver_precision[LEVEL_OUTER] = ms_epsilon;
    MS_solver_precision[LEVEL_PREC]  = prec_ms_epsilon;
    MS_solver_precision[LEVEL_PREC2] = PREC2_MS_EPSILON_FACTOR * prec_ms_epsilon;
    MS_maximum_solver_iterations     = ms_max_iter;

    prec_solver_epsilon   = prec_epsilon;
    prec_solver_max_iter  = prec_max_iter;
    prec2_solver_max_iter = prec_solver_max_iter - PREC2_MAX_ITER_OFFSET;

    /* Populate the partial-fraction tables per level from the shared spectral
       interval. Order 0 stores the linear coefficient c0; order < 0 is inert. */
    for(int lvl = 0; lvl < LEVEL_COUNT; lvl++)
    {
      int order = Zolotarev_order[lvl];
      if(order < 0)
      {
        shifts[lvl] = NULL; numerators[lvl] = NULL; constant_term[lvl] = 0.;
        continue;
      }
      if(order == 0)
      {
        /* order-0 fallback: sign(X) ≈ c0 X, minimax linear on [min,max_eigv].
           Same 'constant_term' role and scaling as the order>=1 tables (it
           multiplies X, with zero poles). */
        shifts[lvl] = NULL; numerators[lvl] = NULL;
        constant_term[lvl] = 2.0 / (ov_params.min_eigv + ov_params.max_eigv);
        if(scaling_factor != 1.0)
          constant_term[lvl] *= 1.0 / sqrt(scaling_factor);
        continue;
      }
      build_zolotarev_tables(lvl, order, scaling_factor);
    }

    /* L1 preconditioning is active if requested (positive tolerance AND a
       positive safety cap) and its order is valid (>= 0, incl. the order-0
       fallback). */
    prec_on = (prec_solver_epsilon > 0) && (prec_solver_max_iter > 0) \
              && (Zolotarev_order[LEVEL_PREC] >= 0);

    /* Second layer is meaningful only if the first layer runs, the L2 order is
       valid (>= 0), and its cap leaves room to iterate. */
    second_layer_on = SECOND_LAYER_REQUESTED
                      && (prec_solver_max_iter > 0)
                      && (Zolotarev_order[LEVEL_PREC2] >= 0)
                      && (prec2_solver_max_iter > 1);

    /* MSCG workspace sized for the larger (outer) Zolotarev order. */
    qpb_mscongrad_init(Zolotarev_order[LEVEL_OUTER]);

  }

  return;
}


void
qpb_overlap_Zolotarev_finalize()
{
  qpb_comm_halo_spinor_field_finalize();
  for(int i=0; i<OVERLAP_NUMB_TEMP_VECS; i++)
    qpb_spinor_field_finalize(ov_temp_vecs[i]);

  for(int i=0; i<MSCG_NUMB_TEMP_VECS; i++)
    qpb_spinor_field_finalize(mscg_temp_vecs[i]);

  if(which_dslash_op == QPB_DSLASH_STANDARD)
    qpb_gauge_field_finalize(*(qpb_gauge_field *)ov_params.gauge_ptr);

  ov_params.initialized = 0;

  qpb_mscongrad_finalize(Zolotarev_order[LEVEL_OUTER]);

  for(int lvl = 0; lvl < LEVEL_COUNT; lvl++)
  {
    if(numerators[lvl]) free(numerators[lvl]);
    if(shifts    [lvl]) free(shifts    [lvl]);
  }

  /* release the deflation subspace */
  if(defl_built)
  {
    for(int i=0; i<defl_k; i++)
      qpb_spinor_field_finalize(defl_V[i]);
    defl_built = 0;
    defl_k = 0;
  }

  return;
}


static void
apply_gamma5_sign(overlap_level_t lvl, qpb_spinor_field y, qpb_spinor_field x)
{
  /* Implements: γ5(sign(X(x))) = γ5(X(c_0 + Sum_{i=1}^{n} c_i/(X^2+σ_i) )),
      with X(x) = γ5(D(x) - ρ*x), at the Zolotarev order of this level.

     With low-mode deflation (defl_built, k > 0):

        γ5 sign(X) x =  D_op( R(X^2) P_high x )         [ high modes, rational ]
                     +  γ5 V sign(W) (V^dag x)          [ low modes, exact     ]

     where P_high x = x - V (V^dag x).  Only the RHS is deflated, so the MSCG
     zeta-recurrence is untouched and every pole keeps its full accuracy. */

  int n = Zolotarev_order[lvl];
  int k = defl_built ? defl_k : 0;

  /* ---- order-0 fallback: sign(X) ≈ c0 X (minimax linear) ----
     Deflated:  γ5 sign(X) x = γ5 V sign(W)(V†x) + c0 (D-ρ) P_high x . */
  if(n == 0)
  {
    qpb_double c0 = constant_term[lvl];

    if(k > 0)
    {
      qpb_spinor_field rhs = ov_temp_vecs[26];
      cdbl a[QPB_DEFL_KMAX];
      for(int i=0; i<k; i++)
        qpb_spinor_xdoty(&a[i], defl_V[i], x);
      qpb_spinor_xeqy(rhs, x);
      for(int i=0; i<k; i++)
        qpb_spinor_axpy(rhs, (qpb_complex){-a[i].re, -a[i].im}, defl_V[i], rhs);

      D_op(y, rhs);                                  /* (D-ρ) P_high x        */
      qpb_spinor_ax(y, (qpb_complex){c0, 0.}, y);    /* c0 (D-ρ) P_high x     */

      cdbl b[QPB_DEFL_KMAX];
      for(int i=0; i<k; i++)
      {
        cdbl s = {0., 0.};
        for(int j=0; j<k; j++)
        {
          cdbl w = defl_signW[i*k+j];
          s.re += w.re*a[j].re - w.im*a[j].im;
          s.im += w.re*a[j].im + w.im*a[j].re;
        }
        b[i] = s;
      }
      qpb_spinor_field low = ov_temp_vecs[0];
      qpb_spinor_field_set_zero(low);
      for(int i=0; i<k; i++)
        qpb_spinor_axpy(low, (qpb_complex){b[i].re, b[i].im}, defl_V[i], low);
      qpb_spinor_gamma5(low, low);
      qpb_spinor_axpy(y, (qpb_complex){1., 0.}, low, y);
    }
    else
    {
      D_op(y, x);
      qpb_spinor_ax(y, (qpb_complex){c0, 0.}, y);    /* c0 (D-ρ) x            */
    }
    return;
  }

  /* ---- order >= 1: deflated Zolotarev sign ---- */

  qpb_spinor_field sum = ov_temp_vecs[0];
  /* MSCG right-hand side: the deflated P_high x when k>0, else x itself. */
  qpb_spinor_field rhs = (k > 0) ? ov_temp_vecs[26] : x;

  qpb_spinor_field yMS[n];
  for(int sigma=0; sigma<n; sigma++)
  {
    yMS[sigma] = mscg_temp_vecs[sigma];
    // It needs to re-initialized to 0 with every call of the function
    qpb_spinor_field_set_zero(yMS[sigma]);
  }

  /* ---- (a) deflate the right-hand side:  a = V^dag x ; rhs = x - V a ---- */
  cdbl a[QPB_DEFL_KMAX];
  if(k > 0)
  {
    for(int i=0; i<k; i++)
      qpb_spinor_xdoty(&a[i], defl_V[i], x);        /* a_i = <V_i, x> */
    qpb_spinor_xeqy(rhs, x);
    for(int i=0; i<k; i++)
      qpb_spinor_axpy(rhs, (qpb_complex){-a[i].re, -a[i].im}, defl_V[i], rhs);
  }

  qpb_double kernel_mass = ov_params.m_bare; // Kernel operator bare mass
  qpb_double kernel_kappa = 1./(2*kernel_mass+8.);

  qpb_mscongrad(yMS, rhs, ov_params.gauge_ptr, ov_params.clover, kernel_kappa, \
    n, shifts[lvl], ov_params.c_sw, MS_solver_precision[lvl], \
    MS_maximum_solver_iterations);

  /* ---- (b) rational part on the (deflated) RHS:  sum = R(X^2) rhs ---- */
  // Initialize sum with the constant term
  qpb_spinor_ax(sum, (qpb_complex) {constant_term[lvl], 0.}, rhs);
  // And then add the rest of the partial fraction terms
  for(int sigma=0; sigma<n; sigma++)
    qpb_spinor_axpy(sum, (qpb_complex) {numerators[lvl][sigma], 0.}, \
                    yMS[sigma], sum);

  D_op(y, sum);   /* high-mode result: γ5 X R(X^2) P_high x */

  /* ---- (c) exact low-mode part:  y += γ5 V sign(W) a ---- */
  if(k > 0)
  {
    cdbl b[QPB_DEFL_KMAX];
    for(int i=0; i<k; i++)
    {
      cdbl s = {0., 0.};
      for(int j=0; j<k; j++)
      {
        cdbl w = defl_signW[i*k+j];
        s.re += w.re*a[j].re - w.im*a[j].im;
        s.im += w.re*a[j].im + w.im*a[j].re;
      }
      b[i] = s;
    }
    /* sum ([0]) is dead after D_op -- reuse it to accumulate the low part */
    qpb_spinor_field low = ov_temp_vecs[0];
    qpb_spinor_field_set_zero(low);
    for(int i=0; i<k; i++)
      qpb_spinor_axpy(low, (qpb_complex){b[i].re, b[i].im}, defl_V[i], low);
    qpb_spinor_gamma5(low, low);
    qpb_spinor_axpy(y, (qpb_complex){1., 0.}, low, y);   /* y += γ5 V sign(W) a */
  }

  return;
}


static void
apply_overlap(overlap_level_t lvl, qpb_spinor_field y, qpb_spinor_field x)
{
  /* Implements:
        Dov,m(x) = (rho+overlap_mass/2)*x + (rho-overlap_mass/2)*g5(sign(X))
     at the Zolotarev order of this level. Preconditioner mass = overlap mass
     (no separate parameter). */

  qpb_spinor_field z = ov_temp_vecs[1];

  qpb_double overlap_mass = ov_params.mass; // Overlap operator Dov,m mass
  qpb_double rho = ov_params.rho;

  qpb_complex a = {rho + 0.5*overlap_mass, 0.};
  qpb_complex b = {rho - 0.5*overlap_mass, 0.};

  apply_gamma5_sign(lvl, z, x);

  qpb_spinor_axpby(y, a, x, b, z);

  return;
}


/* Public single-arg wrappers (callers from main programs). */
void
qpb_gamma5_sign_function_of_X_Zolotarev(qpb_spinor_field y, qpb_spinor_field x)
{
  apply_gamma5_sign(LEVEL_OUTER, y, x);
  return;
}


void
qpb_overlap_Zolotarev(qpb_spinor_field y, qpb_spinor_field x)
{
  apply_overlap(LEVEL_OUTER, y, x);
  return;
}


void
qpb_gamma5_overlap_Zolotarev(qpb_spinor_field y, qpb_spinor_field x)
{
  apply_overlap(LEVEL_OUTER, y, x);
  qpb_spinor_gamma5(y, y);

  return;
}


int
qpb_congrad_overlap_Zolotarev(qpb_spinor_field x, qpb_spinor_field b, \
                                        qpb_double CG_epsilon, int CG_max_iter)
{
  qpb_spinor_field p = ov_temp_vecs[2];
  qpb_spinor_field r = ov_temp_vecs[3];
  qpb_spinor_field y = ov_temp_vecs[4];
  qpb_spinor_field w = ov_temp_vecs[5];
  qpb_spinor_field bprime = ov_temp_vecs[6];

  int n_reeval = 100;
  int n_echo = 100;
  int iters = 0;

  qpb_double res_norm, b_norm;
  qpb_complex_double alpha = {1, 0}, omega = {1, 0};
  qpb_complex_double beta, gamma;

  qpb_spinor_gamma5(w, b);
  qpb_gamma5_overlap_Zolotarev(bprime, w);

  qpb_spinor_xdotx(&b_norm, bprime);

  qpb_spinor_field_set_zero(x);

  /* r0 = bprime - A(x) */
  // qpb_gamma5_overlap_Zolotarev(w, x);
  // qpb_gamma5_overlap_Zolotarev(p, w);
  // qpb_spinor_xmy(r, bprime, p);

  /* Or r0 = bprime for short since x0 = 0 */
  qpb_spinor_xeqy(r, bprime);

  qpb_spinor_xdotx(&gamma.re, r);
  gamma.im = 0;
  res_norm = gamma.re;
  /* p = r0 */
  qpb_spinor_xeqy(p, r);

  qpb_double t = qpb_stop_watch(0);
  for(iters=1; iters<CG_max_iter; iters++)
  {
    // CG stopping criterion
    if(res_norm / b_norm <= CG_epsilon)
    {
      // print("CG stopped at relative residual: %e\n", res_norm / b_norm);
      break;
    }

    /* y = A(p) */
    qpb_gamma5_overlap_Zolotarev(w, p);
    qpb_gamma5_overlap_Zolotarev(y, w);

    /* omega = dot(p, A(p)) */
    qpb_spinor_xdoty(&omega, p, y);

    /* alpha = dot(r, r)/omega */
    alpha = CDEV(gamma, omega);

    /* x <- x + alpha*p */
    qpb_spinor_axpy(x, alpha, p, x);

    if(iters % n_reeval == 0)
    {
      qpb_gamma5_overlap_Zolotarev(w, x);
      qpb_gamma5_overlap_Zolotarev(y, w);
      qpb_spinor_xmy(r, bprime, y);
	  }
    else
	  {
      alpha.re = -CDEVR(gamma, omega);
      alpha.im = -CDEVI(gamma, omega);
      qpb_spinor_axpy(r, alpha, y, r);
	  }
    qpb_spinor_xdotx(&res_norm, r);

    if((iters % n_echo == 0))
	    print(" \t iters = %8d, res = %e\n", iters, res_norm / b_norm);

    beta.re = res_norm / gamma.re;
    beta.im = 0.;
    qpb_spinor_axpy(p, beta, p, r);
    gamma.re = res_norm;
    gamma.im = 0.;
  }

  t = qpb_stop_watch(t);

  qpb_gamma5_overlap_Zolotarev(w, x);
  qpb_gamma5_overlap_Zolotarev(y, w);
  qpb_spinor_xmy(r, bprime, y);
  qpb_spinor_xdotx(&res_norm, r);

  if(iters==CG_max_iter)
  {
    error(" !\n");
    error(" CG *did not* converge, after %d iterations\n", iters);
    error(" residual = %e, relative = %e, t = %g sec\n", res_norm, \
                                                      res_norm / b_norm, t);
    error(" !\n");
    return -1;
  }

  print(" \tAfter %d iters, CG converged, res = %e, relative = %e, "
        "t = %g sec\n",
         iters, res_norm, res_norm / b_norm, t);

  return iters;
}


static int
preconditioner_bicgstab_2(qpb_spinor_field x, qpb_spinor_field b)
{
  /* Level-2 (innermost) preconditioner: right-hand-side solve on D_ov^(n-2),
     applied directly (no further nesting). Iterates until the relative
     recurrence residual ||r||/||b|| <= prec_solver_epsilon, capped at
     prec2_solver_max_iter. Reports iters + recurrence residual (no extra
     operator apply, to keep the dslash accounting clean). */
  qpb_spinor_field r0_hat = ov_temp_vecs[21];
  qpb_spinor_field r      = ov_temp_vecs[22];
  qpb_spinor_field p      = ov_temp_vecs[23];
  qpb_spinor_field u      = ov_temp_vecs[24];
  qpb_spinor_field v      = ov_temp_vecs[25];

  qpb_double res_norm, b_norm;
  qpb_complex_double alpha = {1, 0}, omega = {1, 0};
  qpb_complex_double beta, gamma, rho, zeta;
  int iters;

  qpb_spinor_xdotx(&b_norm, b);
  qpb_spinor_field_set_zero(x);
  qpb_spinor_field_set_zero(p);
  qpb_spinor_field_set_zero(u);

  qpb_spinor_xeqy(r, b);            /* r0 = b - D_ov^prec2·x0 = b */
  qpb_spinor_xeqy(r0_hat, r);

  qpb_spinor_xdotx(&res_norm, r);
  gamma.re = res_norm; gamma.im = 0.;
  rho = gamma;

  for(iters = 1; iters < prec2_solver_max_iter; iters++) {
    if(res_norm / b_norm <= prec_solver_epsilon) break;

    qpb_spinor_xdoty(&gamma, r0_hat, r);
    beta = CMUL(CDEV(gamma, rho), CDEV(alpha, omega));

    omega = CNEGATE(omega);
    qpb_spinor_axpy(p, omega, u, p);      /* p -= omega·u                */
    qpb_spinor_axpy(p, beta, p, r);       /* p  = beta·p + r             */

    apply_overlap(LEVEL_PREC2, u, p);     /* u  = D_ov^prec2 · p         */

    qpb_spinor_xdoty(&beta, r0_hat, u);
    rho   = gamma;
    alpha = CDEV(rho, beta);

    alpha = CNEGATE(alpha);
    qpb_spinor_axpy(r, alpha, u, r);      /* r -= alpha·u   (r ≡ s now)  */

    apply_overlap(LEVEL_PREC2, v, r);     /* v  = D_ov^prec2 · s         */

    qpb_spinor_xdoty(&zeta, v, r);
    qpb_spinor_xdotx(&beta.re, v); beta.im = 0;
    omega = CDEV(zeta, beta);

    alpha = CNEGATE(alpha);
    qpb_spinor_axpy(x, alpha, p, x);      /* x += alpha·p                */
    qpb_spinor_axpy(x, omega, r, x);      /* x += omega·s   (s is in r)  */

    omega = CNEGATE(omega);
    qpb_spinor_axpy(r, omega, v, r);      /* r -= omega·v                */
    omega = CNEGATE(omega);

    qpb_spinor_xdotx(&res_norm, r);       /* recurrence residual         */
  }

  print(" \t\t\tpreconditioner-2 BiCGStab: %d iters, relative residual = %e\n",
        iters, res_norm / b_norm);

  return iters;
}


static void
preconditioner_bicgstab(qpb_spinor_field x, qpb_spinor_field b)
{
  /* Level-1 right-preconditioner solve on D_ov^(n-1): iterates until the
     relative residual ||r||/||b|| <= prec_solver_epsilon, OR until the hard
     safety cap prec_solver_max_iter is reached, whichever comes first. The
     iteration count therefore varies with the right-hand side (non-stationary
     preconditioner; see the note at the prec_solver_* declarations). Its MSCG
     accuracy is MS_solver_precision[LEVEL_PREC].

     When the second layer is enabled, each application of D_ov^(n-1) is itself
     preconditioned by an L2 solve on D_ov^(n-2) (y_pc = K2^{-1}·p, etc.); when
     it is disabled, y_pc ≡ p and z_pc ≡ s, so this reduces exactly to the
     single-layer routine. */
  qpb_spinor_field r0_hat = ov_temp_vecs[14];
  qpb_spinor_field r      = ov_temp_vecs[15];
  qpb_spinor_field p      = ov_temp_vecs[16];
  qpb_spinor_field u      = ov_temp_vecs[17];
  qpb_spinor_field v      = ov_temp_vecs[18];
  qpb_spinor_field y_pc   = ov_temp_vecs[19];   /* K2^{-1}·p */
  qpb_spinor_field z_pc   = ov_temp_vecs[20];   /* K2^{-1}·s */

  qpb_complex_double alpha = {1, 0}, omega = {1, 0};
  qpb_complex_double beta, gamma, rho, zeta;
  qpb_double res_norm, b_norm;
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
    if(res_norm / b_norm <= prec_solver_epsilon) break;

    qpb_spinor_xdoty(&gamma, r0_hat, r);
    beta = CMUL(CDEV(gamma, rho), CDEV(alpha, omega));

    omega = CNEGATE(omega);
    qpb_spinor_axpy(p, omega, u, p);     /* p -= omega·u                 */
    qpb_spinor_axpy(p, beta, p, r);      /* p  = beta·p + r              */

    if(second_layer_on) preconditioner_bicgstab_2(y_pc, p);  /* y_pc = K2⁻¹·p */
    else                qpb_spinor_xeqy(y_pc, p);
    apply_overlap(LEVEL_PREC, u, y_pc);  /* u  = D_ov^(n-1) · y_pc       */

    qpb_spinor_xdoty(&beta, r0_hat, u);
    rho   = gamma;
    alpha = CDEV(rho, beta);

    alpha = CNEGATE(alpha);
    qpb_spinor_axpy(r, alpha, u, r);     /* r -= alpha·u   (r ≡ s now)   */

    if(second_layer_on) preconditioner_bicgstab_2(z_pc, r);  /* z_pc = K2⁻¹·s */
    else                qpb_spinor_xeqy(z_pc, r);
    apply_overlap(LEVEL_PREC, v, z_pc);  /* v  = D_ov^(n-1) · z_pc       */

    qpb_spinor_xdoty(&zeta, v, r);
    qpb_spinor_xdotx(&beta.re, v); beta.im = 0;
    omega = CDEV(zeta, beta);

    alpha = CNEGATE(alpha);
    qpb_spinor_axpy(x, alpha, y_pc, x);  /* x += alpha·y_pc  (NOT alpha·p) */
    qpb_spinor_axpy(x, omega, z_pc, x);  /* x += omega·z_pc  (NOT omega·s) */

    omega = CNEGATE(omega);
    qpb_spinor_axpy(r, omega, v, r);     /* r -= omega·v                 */
    omega = CNEGATE(omega);

    qpb_spinor_xdotx(&res_norm, r);      /* recurrence residual for exit test */
  }

  /* Explicit final (true) residual, reported for diagnostics — the
     recurrence residual used in the exit test can drift from b - D_ov^prec·x. */
  apply_overlap(LEVEL_PREC, u, x);
  qpb_spinor_xmy(r, b, u);
  qpb_spinor_xdotx(&res_norm, r);
  print(" \t\tpreconditioner BiCGStab: %d iters, relative residual = %e\n",
        iters, res_norm / b_norm);

  return;
}


int
qpb_bicgstab_overlap_Zolotarev(qpb_spinor_field x, qpb_spinor_field b, \
                               qpb_double epsilon, int max_iter)
{
  qpb_spinor_field r0_hat = ov_temp_vecs[7];
  qpb_spinor_field r      = ov_temp_vecs[8];
  qpb_spinor_field p      = ov_temp_vecs[9];
  qpb_spinor_field u      = ov_temp_vecs[10];
  qpb_spinor_field v      = ov_temp_vecs[11];
  qpb_spinor_field y_pc   = ov_temp_vecs[12];   /* K^{-1}·p */
  qpb_spinor_field z_pc   = ov_temp_vecs[13];   /* K^{-1}·s */

  int n_reeval = 100, n_echo = 100, iters = 0;
  qpb_double res_norm, b_norm;
  qpb_complex_double alpha = {1, 0}, omega = {1, 0};
  qpb_complex_double beta, gamma, rho, zeta;

  /* Select the preconditioner path once, outside the loop. */
  enum { PREC_NONE, PREC_BICGSTAB } prec_path;
  prec_path = prec_on ? PREC_BICGSTAB : PREC_NONE;

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
    case PREC_NONE:     qpb_spinor_xeqy(y_pc, p);              break;
    case PREC_BICGSTAB: preconditioner_bicgstab(y_pc, p);      break;
    }

    apply_overlap(LEVEL_OUTER, u, y_pc);          /* u = D_ov · y_pc */

    qpb_spinor_xdoty(&beta, r0_hat, u);
    rho   = gamma;
    alpha = CDEV(rho, beta);

    alpha = CNEGATE(alpha);
    qpb_spinor_axpy(r, alpha, u, r);              /* r ≡ s now */

    /* z_pc = K^{-1}·s */
    switch(prec_path) {
    case PREC_NONE:     qpb_spinor_xeqy(z_pc, r);              break;
    case PREC_BICGSTAB: preconditioner_bicgstab(z_pc, r);      break;
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
