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

#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_sf_ellint.h>
#include <gsl/gsl_sf_elljac.h>


#define OVERLAP_NUMB_TEMP_VECS 7
#define MSCG_NUMB_TEMP_VECS 20


static qpb_spinor_field ov_temp_vecs[OVERLAP_NUMB_TEMP_VECS];
static qpb_spinor_field mscg_temp_vecs[MSCG_NUMB_TEMP_VECS];

static qpb_overlap_params ov_params;

static int Zolotarev_order;
static qpb_double MS_solver_precision;
static int MS_maximum_solver_iterations;

static qpb_double *numerators;
static qpb_double *shifts;
static qpb_double constant_term;


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
  algorithm. */

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

/* ---------------- RAYLEIGH QUOTIENT ITERATION REFINEMENT ---------------- */

/* The smallest Ritz value returned by qpb_extreme_eigenvalues_of_X_squared is
a variational *upper* bound on the smallest eigenvalue of X^2: Lanczos
approaches it from above and can still be loose at the iteration counts used
here. The routines below refine that estimate with a warm start by plain
inverse iteration, followed by Rayleigh quotient iteration, which is cubically
convergent once the vector is aligned with the target eigenvector. */

/* Tunables. These are deliberately kept here rather than in the input file:
the point of the refinement is its improvement-vs-cost curve, and these are
the knobs that curve is traced with. */
static const int RQI_n_warm = 3;   /* inverse iteration steps (warm start)  */
static const int RQI_n_rqi = 10;   /* Rayleigh quotient iteration steps     */

/* Stage 2 shifts *below* the current estimate rather than exactly at it:
sigma = (1 - RQI_shift_backoff)*lambda^2. Shifting exactly at lambda^2 makes
p^dag A p = RQ(v) - sigma vanish identically at CG's first iteration, an
exact algebraic cancellation rather than roundoff, which is what previously
killed every shifted solve after at most one step. Backing off keeps
X^2 - sigma positive definite whenever the backoff exceeds the current
relative error, so plain CG is rigorously valid; convergence degrades from
cubic to fast-linear, which is still far better than stage 1 and, unlike
cubic-in-principle, actually survives more than one step. The pAp guard
below stays as a backstop for the steps where the error still exceeds the
backoff. */
static const qpb_double RQI_shift_backoff = 1e-2;

/* Solver tolerances follow the qpb convention throughout: 'epsilon' is
compared against the *squared* relative residual, so 1e-3 below is a true
relative residual of ~3e-2. Stage 1 now only has to get the vector aligned
enough for stage 2 to take over, not converged, since stage 2 can take
several working steps. Stage 2's operator is near-singular by construction
and is meant to be solved loosely, with the iteration cap rather than the
tolerance doing the stopping. */
static const qpb_double RQI_warm_epsilon = 1e-3;
static const int RQI_warm_max_iters = 1000;
static const qpb_double RQI_rqi_epsilon = 1e-4;
static const int RQI_rqi_max_iters = 100;

/* Stop when the Rayleigh quotient stops moving. Its accuracy is limited by
rounding in the matrix-vector product to about eps_machine*||X^2||/lambda^2,
i.e. ~1e-12 relative here, so this threshold sits just above the noise floor. */
static const qpb_double RQI_lambda2_tolerance = 1e-11;

/* Guard on p^dag A p in the shifted solve. (X^2 - lambda^2) always carries one
non-positive eigenvalue, since the Rayleigh quotient is an upper bound on
lambda_min^2 for *every* vector, so the operator is indefinite from the first
Rayleigh quotient iteration step and not merely near convergence. CG can then
either flip the sign of its step length (p^dag A p < 0) or break down outright
(p^dag A p ~ 0); the single test below catches both. */
static const qpb_double RQI_pAp_floor = 1e-14;


INLINE void
X_op(qpb_spinor_field y, qpb_spinor_field x)
{
  /* Implements X = g5*(D - rho). Note that ov_params.gauge_ptr already
  carries the time boundary condition, applied once in
  qpb_overlap_Zolotarev_init, so unlike qpb_congrad() nothing is done to the
  gauge field here. */

  void *dslash_args[4];

  dslash_args[0] = ov_params.gauge_ptr;
  dslash_args[1] = &ov_params.m_bare;
  dslash_args[2] = &ov_params.clover;
  dslash_args[3] = &ov_params.c_sw;

  ov_params.g5_dslash_op(y, x, dslash_args);

  return;
}


static int
normalize_spinor(qpb_spinor_field v, qpb_spinor_field x)
{
  /* v <- x/||x||. Returns 0 if x has no usable norm, which happens if a
  solve bailed out before it ever updated its solution vector. The test is
  written so that a NaN norm fails it too. */

  qpb_double norm;
  qpb_spinor_xdotx(&norm, x);

  if(!(norm > 0.))
    return 0;

  qpb_spinor_ax(v, (qpb_complex){1./sqrt(norm), 0.}, x);

  return 1;
}


static int
RQI_congrad(qpb_spinor_field x, qpb_spinor_field b, qpb_double shift, \
            qpb_double epsilon, int max_iters, \
            qpb_spinor_field p, qpb_spinor_field r, qpb_spinor_field y, \
            qpb_spinor_field w, int *bailed)
{
  /* CG on (X^2 + shift), starting from x = 0. Returns the number of
  iterations performed and sets '*bailed' if the p^dag A p guard stopped it.
  Reaching 'max_iters' is not treated as an error: for the shifted systems of
  stage 2 that is the expected exit, and the current iterate is still a
  perfectly good direction. */

  qpb_double res_norm, b_norm, p_norm;
  qpb_complex_double alpha, omega, gamma, beta;
  qpb_complex c_shift = (qpb_complex){shift, 0.};

  *bailed = 0;

  /* x = 0, and therefore r = b - (X^2 + shift) x = b */
  qpb_spinor_field_set_zero(x);
  qpb_spinor_xeqy(r, b);
  qpb_spinor_xeqy(p, b);

  qpb_spinor_xdotx(&b_norm, b);
  res_norm = b_norm;
  gamma = (qpb_complex_double){b_norm, 0.};

  int iters;
  for(iters=1; iters<max_iters; iters++)
  {
    if(res_norm / b_norm <= epsilon)
      break;

    /* w = (X^2 + shift) p */
    X_op(y, p);
    X_op(w, y);
    qpb_spinor_axpy(w, c_shift, p, w);

    qpb_spinor_xdoty(&omega, p, w);
    qpb_spinor_xdotx(&p_norm, p);

    if(omega.re <= RQI_pAp_floor*p_norm)
    {
      print("   RQI: CG step %d bailed on p^dag A p = %e\n", iters, omega.re);
      *bailed = 1;
      break;
    }

    alpha = CDEV(gamma, omega);
    qpb_spinor_axpy(x, alpha, p, x);

    alpha.re = -CDEVR(gamma, omega);
    alpha.im = -CDEVI(gamma, omega);
    qpb_spinor_axpy(r, alpha, w, r);

    qpb_spinor_xdotx(&res_norm, r);

    beta.re = res_norm / gamma.re;
    beta.im = 0.;
    qpb_spinor_axpy(p, beta, p, r);
    gamma.re = res_norm;
    gamma.im = 0.;
  }

  return iters;
}


static qpb_double
refine_min_eigenvalue_RQI(qpb_double lambda2_in, int *n_warm_solves, \
                          int *n_rqi_steps, int *n_cg_iters)
{
  /* Refine the smallest eigenvalue of X^2 by inverse iteration followed by
  Rayleigh quotient iteration.
    lambda2_in : Ritz estimate from qpb_extreme_eigenvalues_of_X_squared
    returns    : refined lambda_min^2, never larger than lambda2_in
  The cost is reported through the three counters. */

  qpb_spinor_field p = qpb_spinor_field_init();
  qpb_spinor_field r = qpb_spinor_field_init();
  qpb_spinor_field y = qpb_spinor_field_init();
  qpb_spinor_field w = qpb_spinor_field_init();
  qpb_spinor_field v = qpb_spinor_field_init();
  qpb_spinor_field z = qpb_spinor_field_init();
  qpb_spinor_field t = qpb_spinor_field_init();

  int bailed;

  *n_warm_solves = 0;
  *n_rqi_steps = 0;
  *n_cg_iters = 0;

  /* ------- stage 1: warm start by plain inverse iteration on X^2 -------
  Unshifted, and deliberately so: applying (X^2)^{-1} already amplifies the
  component belonging to the smallest eigenvalue, since that component is
  divided by the smallest number, so no prior knowledge of lambda_min^2 is
  needed. It also keeps the operator SPD, which shifting by the Ritz value
  would not: Ritz values approach from above, so X^2 - lambda_Ritz^2 is
  typically indefinite. Convergence here is only linear, but this stage only
  has to produce a vector reasonably aligned with the target eigenvector. */

  qpb_spinor_field_set_random(v);
  if(!normalize_spinor(v, v))
    {
      /* Cannot happen in practice; the tracker below still guarantees that
      the routine returns no worse than the Ritz value it was given. */
      error(" RQI: random start vector has zero norm\n");
    }

  for(int k=0; k<RQI_n_warm; k++)
  {
    *n_cg_iters += RQI_congrad(z, v, 0., RQI_warm_epsilon, \
                               RQI_warm_max_iters, p, r, y, w, &bailed);
    *n_warm_solves += 1;

    if(!normalize_spinor(v, z))
    {
      print("   RQI: warm-start solve %d produced no usable vector,"
            " stopping the warm start\n", k+1);
      break;
    }
  }

  /* Stage 1 runs on the unshifted operator, so its own Rayleigh quotient
  carries no indefiniteness risk to evaluate and is directly comparable to
  the Ritz value. In practice it is often *better* than a loosely-converged
  Lanczos estimate: 10 solves at a decent tolerance can outperform a Lanczos
  run cut short after a similar number of iterations. Computing it costs one
  more matrix-vector pair on vectors already in hand, and folding it into the
  running best means a weak Ritz value can never mask a gain stage 1 already
  paid for. */
  qpb_complex_double warm_numerator;
  qpb_double warm_denominator;

  X_op(y, v);
  X_op(t, y);
  qpb_spinor_xdoty(&warm_numerator, v, t);
  qpb_spinor_xdotx(&warm_denominator, v);

  qpb_double lambda2_warm = warm_numerator.re / warm_denominator;

  print("   RQI: warm start  lambda^2 = %.16e  (Ritz was %.16e)\n", \
        lambda2_warm, lambda2_in);

  /* ---------------- stage 2: Rayleigh quotient iteration ----------------
  The shift is re-derived from the current vector at every step, but backed
  off below it rather than placed exactly on it: shifting exactly at the
  Rayleigh quotient is what textbook RQI does to earn cubic convergence, and
  it is also precisely what makes the first CG search direction satisfy
  p^dag A p = 0 identically, so plain CG cannot take a single step. Backing
  off trades the cubic rate for a fast-linear one that CG can actually
  deliver. */

  /* Every Rayleigh quotient is an upper bound on lambda_min^2, so the
  smallest one seen is the best estimate available. Seeding this with the
  Ritz value means the returned number can never be worse than the one the
  refinement started from, whatever the shifted solves do; folding in the
  warm-start value below extends that guarantee to stage 1 as well. */
  qpb_double lambda2_best = lambda2_in;
  if(lambda2_warm < lambda2_best)
    lambda2_best = lambda2_warm;

  /* Start stage 2 from whichever of the two is already tighter, so the
  first shifted solve is never handed a shift worse than what stage 1 found
  for free: that gap is exactly what previously showed up as an immediate,
  uninformative p^dag A p < 0 bailout. */
  qpb_double lambda2 = lambda2_best;

  print("   RQI: step  0  lambda^2 = %.16e  (starting shift)\n", lambda2);

  for(int i=0; i<RQI_n_rqi; i++)
  {
    qpb_double sigma = (1. - RQI_shift_backoff)*lambda2;

    int iters = RQI_congrad(z, v, -sigma, RQI_rqi_epsilon, \
                            RQI_rqi_max_iters, p, r, y, w, &bailed);
    *n_cg_iters += iters;
    *n_rqi_steps += 1;

    /* How the solve ended is the diagnostic that says whether the backoff
    was generous enough: 'cap' and 'converged' both mean CG stayed valid
    throughout, 'bailed' means the guard still had to catch it. */
    const char *exit_reason = bailed ? "bailed" : \
                              (iters >= RQI_rqi_max_iters ? "cap" : "converged");

    if(!normalize_spinor(v, z))
    {
      print("   RQI: step %2d produced no usable vector (%s), stopping\n", \
            i+1, exit_reason);
      break;
    }

    /* Rayleigh quotient lambda^2 = (v, X^2 v)/(v, v) */
    qpb_complex_double numerator;
    qpb_double denominator;

    X_op(y, v);
    X_op(t, y);
    qpb_spinor_xdoty(&numerator, v, t);
    qpb_spinor_xdotx(&denominator, v);

    qpb_double lambda2_new = numerator.re / denominator;

    print("   RQI: step %2d  lambda^2 = %.16e  "
          "(sigma = %.10e, CG iters = %4d, %s)\n", \
          i+1, lambda2_new, sigma, iters, exit_reason);

    if(lambda2_new < lambda2_best)
      lambda2_best = lambda2_new;

    qpb_double change = fabs(lambda2_new - lambda2) / fabs(lambda2_new);
    lambda2 = lambda2_new;

    if(change < RQI_lambda2_tolerance)
      break;
  }

  qpb_spinor_field_finalize(p);
  qpb_spinor_field_finalize(r);
  qpb_spinor_field_finalize(y);
  qpb_spinor_field_finalize(w);
  qpb_spinor_field_finalize(v);
  qpb_spinor_field_finalize(z);
  qpb_spinor_field_finalize(t);

  return lambda2_best;
}


/* ------------------------ MATRIX-VECTOR FUNCTIONS ------------------------ */

void
qpb_overlap_Zolotarev_init(void * gauge, qpb_clover_term clover, \
          int Zol_iters, qpb_double rho, \
          qpb_double c_sw, qpb_double mass, qpb_double scaling_factor, \
          qpb_double ms_epsilon, int ms_max_iter, \
          qpb_double Lanczos_epsilon, int Lanczos_max_iters, \
          qpb_double delta_max, qpb_double delta_min, \
          int RQI_refinement_mode)
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

    /* --------------------- extreme eigenvalues of X^2 --------------------- */

    qpb_double min_eigv_squared;
    qpb_double max_eigv_squared;

    /* First the the extrema of the eigenvalues spectrum of X^2,
    X = g5*(D - rho), are calculated and are stored inside the
    'min_eigv_squared' and 'max_eigv_squared'variables correspondingly. */
    int Lanczos_iters = qpb_extreme_eigenvalues_of_X_squared(&min_eigv_squared,\
                      &max_eigv_squared, Lanczos_epsilon, Lanczos_max_iters);
    print(" Total number of Lanczos algorithm iterations = %d\n", \
                                                                Lanczos_iters);

    /* The Ritz estimate of lambda_min^2 is refined here, before any of the
    downstream modifications are applied, so that what is reported is directly
    comparable to the Ritz value it replaces. */
    if(RQI_refinement_mode != QPB_RQI_OFF)
    {
      int n_warm_solves, n_rqi_steps, n_cg_iters;

      qpb_double t_RQI = qpb_stop_watch(0);
      qpb_double refined_min_eigv_squared = refine_min_eigenvalue_RQI(\
              min_eigv_squared, &n_warm_solves, &n_rqi_steps, &n_cg_iters);
      t_RQI = qpb_stop_watch(t_RQI);

      print(" lambda_min^2 (Ritz)                      = %.10e\n", \
                                                            min_eigv_squared);
      print(" lambda_min^2 (after RQI)                 = %.10e\n", \
                                                    refined_min_eigv_squared);
      print(" relative change                          = %.2e\n", \
              fabs(refined_min_eigv_squared - min_eigv_squared)\
                                              / fabs(min_eigv_squared));
      print(" RQI: warm-start solves / RQI steps       = %d / %d\n", \
                                                n_warm_solves, n_rqi_steps);
      print(" RQI: total CG iterations                 = %d\n", n_cg_iters);
      print(" RQI: Lanczos iterations, for comparison  = %d\n", Lanczos_iters);
      print(" RQI: time                                = %g secs\n", t_RQI);

      if(RQI_refinement_mode == QPB_RQI_APPLY)
      {
        min_eigv_squared = refined_min_eigv_squared;
      }
      else
      {
        print(" RQI: measure-only, the Ritz value is used downstream\n");
      }
    }

    /* If requested the extreme eigenvalues are modified accordingly */
    if (delta_min != 1.0)
      min_eigv_squared *= delta_min;
    if (delta_max != 1.0)
      max_eigv_squared *= delta_max;
    
    print(" Min eigenvalue squared = %.16f\n", min_eigv_squared);
    print(" Max eigenvalue squared = %.16f\n", max_eigv_squared);

    /* And then their square root value is stored inside the 'min_eigv' and
    'max_eigv' attributes of the 'ov_params' struct. */
    
    ov_params.min_eigv = sqrt(min_eigv_squared);
    ov_params.max_eigv = sqrt(max_eigv_squared);

    /* ----------------------- expansion coefficients ----------------------- */

    Zolotarev_order = Zol_iters;
    MS_solver_precision = ms_epsilon;
    MS_maximum_solver_iterations = ms_max_iter;

    /* Calculate the numerical terms of the partial fraction expansion */
    shifts = qpb_alloc(sizeof(qpb_double)*Zolotarev_order);
    numerators = qpb_alloc(sizeof(qpb_double)*Zolotarev_order);

    /* Compute elliptic parameters */
    qpb_double k_squared = min_eigv_squared/max_eigv_squared;
    qpb_double k_prime = sqrt(1.0 - k_squared);
    
    /* Compute complete elliptic integrals */
    qpb_double K_prime = gsl_sf_ellint_Kcomp(k_prime, GSL_PREC_DOUBLE);

    /* Allocate arrays */
    qpb_double *c = qpb_alloc(sizeof(qpb_double)*2*Zolotarev_order);
    qpb_double *b = qpb_alloc(sizeof(qpb_double)*Zolotarev_order);
    
    /* Compute coefficients */
    calculate_c_coefficients(c, Zolotarev_order, k_prime, K_prime);
    calculate_b_coefficients(b, c, Zolotarev_order);
    
    /* Compute normalization constant */
    qpb_double normalization_constant = calculate_normalization_constant(
                                  c, Zolotarev_order, \
                                  ov_params.min_eigv, ov_params.max_eigv);

    constant_term = normalization_constant / ov_params.min_eigv;

    for(int i=0; i<Zolotarev_order; i++)
    {
      shifts[i] = c[2*i] * ov_params.min_eigv * ov_params.min_eigv;
      numerators[i] = normalization_constant * b[i] * ov_params.min_eigv;
      // print("numerator[%d] = %.25f, shift[%d] = %.25f\n", i, numerators[i], \
                                                              i, shifts[i]);
    }

    // Modify the numerical constants of the partial fraction expansions using
    // the scaling parameter
    if (scaling_factor != 1.0)
    {
      constant_term *= 1/sqrt(scaling_factor);
      for(int i=0; i<Zolotarev_order; i++)
      {
        numerators[i] *= sqrt(scaling_factor);
        shifts[i] *= scaling_factor;
      }
    }

    qpb_mscongrad_init(Zolotarev_order);

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
  
  qpb_mscongrad_finalize(Zolotarev_order);

  free(numerators);
  free(shifts);
  
  return;
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


void
qpb_gamma5_sign_function_of_X_Zolotarev(qpb_spinor_field y, qpb_spinor_field x)
{
  /* Implements: γ5(sign(X(x))) = γ5(X(c_0 + Sum_{i=1}^{n} c_i/(X^2+σ_i) )),
      with X(x) = γ5(D(x) - ρ*x) . */

  qpb_spinor_field sum = ov_temp_vecs[0];

  qpb_spinor_field yMS[Zolotarev_order];
  for(int sigma=0; sigma<Zolotarev_order; sigma++)
  {
    yMS[sigma] = mscg_temp_vecs[sigma];
    // It needs to re-initialized to 0 with every call of the function
    qpb_spinor_field_set_zero(yMS[sigma]);
  }

  qpb_double kernel_mass = ov_params.m_bare; // Kernel operator bare mass
  qpb_double kernel_kappa = 1./(2*kernel_mass+8.);

  qpb_mscongrad(yMS, x, ov_params.gauge_ptr, ov_params.clover, kernel_kappa, \
    Zolotarev_order, shifts, ov_params.c_sw, MS_solver_precision, \
    MS_maximum_solver_iterations);

  // Initialize sum with the constant term
  qpb_spinor_ax(sum, (qpb_complex) {constant_term, 0.}, x);
  // And then add the rest of the partial fraction terms
  for(int sigma=0; sigma<Zolotarev_order; sigma++)
    qpb_spinor_axpy(sum, (qpb_complex) {numerators[sigma], 0.}, yMS[sigma], sum);

  D_op(y, sum);

  return;
}


void
qpb_overlap_Zolotarev(qpb_spinor_field y, qpb_spinor_field x)
{
  /* Implements:
        Dov,m(x) = (rho+overlap_mass/2)*x+ (rho-overlap_mass/2)*g5(sign(X))
  */
  
  qpb_spinor_field z = ov_temp_vecs[1];

  qpb_double overlap_mass = ov_params.mass; // Overlap operator Dov,m mass
  qpb_double rho = ov_params.rho;

  qpb_complex a = {rho + 0.5*overlap_mass, 0.};
  qpb_complex b = {rho - 0.5*overlap_mass, 0.};

  qpb_gamma5_sign_function_of_X_Zolotarev(z, x);

  qpb_spinor_axpby(y, a, x, b, z);

  return;
}


void
qpb_gamma5_overlap_Zolotarev(qpb_spinor_field y, qpb_spinor_field x)
{
  qpb_overlap_Zolotarev(y, x);
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
