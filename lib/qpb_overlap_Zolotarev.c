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

        /* sn2/(1-sn2) cancels catastrophically as sn -> 1 (large beta/alpha,
        e.g. relative error ~0.26 at b ~ 1e7); use sn^2+cn^2=1 instead. */
        c[i-1] = (sn * sn) / (cn * cn);
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


/* Build the (c, b, norm) triple of a degree-(n,n) Zolotarev approximation
on |X| in [min_eigv, max_eigv]. Factored out so it can be called repeatedly
by select_delta_min() (once per trial alpha) as well as once by
qpb_overlap_Zolotarev_init(), instead of duplicating the construction. */
static void
build_Zolotarev_coeffs(qpb_double min_eigv, qpb_double max_eigv, int n,
                       qpb_double *c, qpb_double *b, qpb_double *norm)
{
  qpb_double k_squared = (min_eigv*min_eigv)/(max_eigv*max_eigv);
  qpb_double k_prime   = sqrt(1.0 - k_squared);
  qpb_double K_prime   = gsl_sf_ellint_Kcomp(k_prime, GSL_PREC_DOUBLE);
  calculate_c_coefficients(c, n, k_prime, K_prime);
  calculate_b_coefficients(b, c, n);
  *norm = calculate_normalization_constant(c, n, min_eigv, max_eigv);
}


/* Z_n at a PHYSICAL argument x, for the coefficients (c, norm) built by
build_Zolotarev_coeffs() with the given min_eigv: sign_function_product_form
expects the rescaled argument x/min_eigv, the division is internal here so
callers always pass physical spectral points. */
static qpb_double
Zolotarev_at(qpb_double x, qpb_double min_eigv, int n,
             qpb_double *c, qpb_double norm)
{
  return norm * sign_function_product_form(x / min_eigv, n, c);
}


/* Select the largest delta_min = alpha^2/min_eigv_sq_raw, capped at 0.5, such
that the resulting Zolotarev approximation keeps a healthy margin
sigma_min = C - R*Z_n(x) >= theta_window*am across the window
[rho - 1.1*lamR, rho - 0.9*lamR], where lamR is the leftmost real kernel mode
(see ZOLOTAREV_DELTA_MIN_SPEC.md, ZOLOTAREV_DELTA_MIN_PATCH.md and
ZOLOTAREV_PATCH_3.md for the physics). Returns -1.0 if no trial alpha in the
scan satisfies the margin everywhere in the window: that means no lower
bound is simultaneously safe and accurate at this order, and the caller
must raise n rather than fall back to a default.

This window is intentionally narrow, testing only the immediate
neighbourhood of the resonance point rho - delta: widening it to the full
[lambda_min, rho - delta] range makes the interior 20-point scan also catch
genuine mid-range equioscillation overshoot at low n (e.g. n=2), which can
reject alpha values that are otherwise fine at rho - delta itself - at
n=2, am=0.05 on the reference config this rejects every trial alpha, where
the narrow window here finds sigma_min/am ~ 1.0. That broader overshoot is
real and worth surfacing, but as a report, not a rejection criterion -
see print_overshoot_regions(), which scans the wide range separately.

theta_window is named to disambiguate it from the unrelated theta_gate used
by the d_Z quality gate in qpb_overlap_Zolotarev_init: this one controls how
aggressively alpha is reduced when the danger window is compromised, and
stays at the decided value 0.5; theta_gate controls a chiral-quality report
on the final choice and is a different constant entirely. Do not conflate
the two.

The margin test is signed (C - R*Z, not C - R*|Z|), so it still rejects the
index-flipped case (Z > C/R gives a negative margin) without needing a
separate Z < 1 test: the old Z >= 1.0 check was stricter than the physics
requires and rejected healthy Z in (1, C/R] for no reason, driving alpha
down by orders of magnitude at moderate n. Z_n is allowed to exceed 1 here
(the equioscillation overshoot is normal for this normalization); the
caller reports where it does via the overshoot-region diagnostic instead.

The delta_min <= 0.5 cap preserves the *original*, resonance-independent
reason for shrinking below the raw Lanczos estimate: Lanczos converges to
lambda_min from above, so the true spectrum can extend below alpha, and 0.5
is the conventional safety margin for that. Without the cap, the tightened
margin test above would return exactly 1.0 (no reduction at all) whenever
the resonance itself isn't a problem at the requested n.

Do not replace the full ascending scan with a bisection or a "first pass,
scanning downward" search: for n >= 2, Z_n(rho - delta) is not monotonic in
alpha (multiple crossings), so only scanning the whole range and keeping the
largest passing alpha is robust to arbitrary lobe structure. */
static qpb_double
select_delta_min(qpb_double min_eigv_sq_raw, qpb_double max_eigv,
                 qpb_double rho, qpb_double am, qpb_double lamR,
                 int n, qpb_double theta_window)
{
  qpb_double lmin = sqrt(min_eigv_sq_raw);
  qpb_double alpha_max = lmin * sqrt(0.5);        /* delta_min <= 0.5 */
  qpb_double C = rho + 0.5*am, R = rho - 0.5*am;
  qpb_double xlo = rho - 1.1*lamR, xhi = rho - 0.9*lamR;
  qpb_double *c = qpb_alloc(sizeof(qpb_double)*2*n);
  qpb_double *b = qpb_alloc(sizeof(qpb_double)*n);
  qpb_double best = -1.0, norm;

  for(int j=0; j<200; j++)                 /* ascending in alpha */
  {
    qpb_double alpha = alpha_max * pow(10.0, -4.0 + 4.0*j/199.0);
    build_Zolotarev_coeffs(alpha, max_eigv, n, c, b, &norm);

    int ok = 1;
    for(int i=0; i<20 && ok; i++)
    {
      qpb_double x = xlo + (xhi - xlo)*i/19.0;
      qpb_double Z = Zolotarev_at(x, alpha, n, c, norm);
      if (C - R*Z < theta_window*am)   ok = 0;   /* signed margin: also    */
    }                                             /* rejects index flips   */
    if (ok) best = alpha;                  /* keep the LARGEST passing alpha */
  }

  free(c); free(b);
  return (best > 0.0) ? (best*best)/min_eigv_sq_raw : -1.0;
}


/* Report the sub-intervals of [xlo, xhi] where Z_n(x) > C/R, i.e. where a
real kernel mode would have a flipped index. Called with the WIDE range
[lambda_min, rho - delta], deliberately wider than select_delta_min's own
narrow acceptance window [rho - 1.1*delta, rho - 0.9*delta]: the selector
only needs to guard the immediate neighbourhood of the resonance point, but
a genuine overshoot further out (as seen at n=2) is real and still worth
surfacing even though it isn't grounds for rejecting alpha. Diagnostic-only,
always printed regardless of whether delta_min needed reducing (see
ZOLOTAREV_DELTA_MIN_PATCH.md #6 and ZOLOTAREV_PATCH_3.md). */
static void
print_overshoot_regions(qpb_double xlo, qpb_double xhi, qpb_double alpha,
                        qpb_double C, qpb_double R, int n, qpb_double *c,
                        qpb_double norm)
{
  const int npts = 500;
  qpb_double threshold = C/R;
  int printed_any = 0;
  qpb_double region_lo = 0.0;
  int in_region = 0;

  for(int i=0; i<npts; i++)
  {
    qpb_double x = xlo + (xhi - xlo)*(i+1)/npts;
    qpb_double Z = Zolotarev_at(x, alpha, n, c, norm);
    int over = (Z > threshold);

    if(over && !in_region)
    {
      region_lo = x;
      in_region = 1;
    }
    else if(!over && in_region)
    {
      print(" overshoot region (Z_n > C/R) on [lambda_min, rho-delta] : "
            "[%.6f, %.6f]\n", region_lo, x);
      printed_any = 1;
      in_region = 0;
    }
  }
  if(in_region)
  {
    print(" overshoot region (Z_n > C/R) on [lambda_min, rho-delta] : "
          "[%.6f, %.6f]\n", region_lo, xhi);
    printed_any = 1;
  }
  if(!printed_any)
    print(" overshoot region (Z_n > C/R) on [lambda_min, rho-delta] : none\n");
}


/* Equioscillation amplitude d_Z = max_{[alpha,beta]} (Z_n(x) - 1). This is
independent of the window/margin test above: that test only guards
invertibility and index correctness near rho - delta, while chiral symmetry
breaking is governed by d_Z over the *whole* interval. The two are
demonstrably different quantities - a run can pass the window test with
sigma_min/am ~ 1 while still carrying a large residual mass from d_Z (see
ZOLOTAREV_PATCH_3.md). Report-only: called once on the final selected alpha,
not inside the delta_min scan. */
static qpb_double
equioscillation_amplitude(qpb_double alpha, qpb_double beta, int n)
{
  qpb_double *c = qpb_alloc(sizeof(qpb_double)*2*n);
  qpb_double *b = qpb_alloc(sizeof(qpb_double)*n);
  qpb_double norm, d_Z = -1.0;

  build_Zolotarev_coeffs(alpha, beta, n, c, b, &norm);
  for(int j=0; j<200; j++)
  {
    qpb_double x = alpha * pow(beta/alpha, (qpb_double)j/199.0);
    qpb_double dev = Zolotarev_at(x, alpha, n, c, norm) - 1.0;
    if (dev > d_Z) d_Z = dev;
  }
  free(c); free(b);
  return d_Z;
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
qpb_extreme_eigenvalues_of_X_squared(qpb_double *min_eigv,
    qpb_double *max_eigv, qpb_double bare_mass,
    qpb_double Lanczos_epsilon, int max_iters, int min_iters)
{
  /* It calculates the extreme eigenvalues of the eigenvalue spectrum
  of H^2, H ≡ γ5*Kernel(x), with: Kernel(x) = (a*D - ρ)(x), using the Lanczos
  algorithm. 'min_iters' is a floor on the number of iterations before the
  convergence check is allowed to stop the loop: Lanczos approaches
  lambda_min from above and can sit on a plateau, so a loose tolerance
  combined with a small iteration count can stop early on a not-yet-converged
  estimate (pass 0 for the historical no-floor behavior). */

  qpb_lanczos_init();

  qpb_clover_term clover_term = ov_params.clover;
  qpb_double c_sw = ov_params.c_sw;
  qpb_double mass = bare_mass;
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
    if((i+1 >= min_iters) && dlambda < Lanczos_epsilon*0.5)
      break;
    lambda0 = lambda;
  }

  *min_eigv = (qpb_double) eig[0];
  *max_eigv = (qpb_double) eig[i-1];

  free(a);
  free(b);
  free(eig);

  qpb_lanczos_finalize();

  return i;
}

/* ------------------------ MATRIX-VECTOR FUNCTIONS ------------------------ */

void
qpb_overlap_Zolotarev_init(void * gauge, qpb_clover_term clover, \
          int Zol_iters, qpb_double rho, \
          qpb_double c_sw, qpb_double mass, qpb_double scaling_factor, \
          qpb_double ms_epsilon, int ms_max_iter, \
          qpb_double Lanczos_epsilon, int Lanczos_max_iters, \
          qpb_double delta_max, qpb_double delta_min, \
          qpb_double kernel_delta)
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
                      &max_eigv_squared, ov_params.m_bare,  // Kernel operator mass set at -rho
                      Lanczos_epsilon, Lanczos_max_iters, 0);
    print(" Total number of Lanczos algorithm iterations = %d\n", \
                                                                Lanczos_iters);
    /* delta_min is kept only for source compatibility with existing callers;
    a fixed multiplier here can drive the Zolotarev approximation into
    resonance with the leftmost kernel mode (see
    ZOLOTAREV_DELTA_MIN_SPEC.md), so it is no longer applied. The safe
    replacement, 'selected_delta_min' below, is computed automatically. */
    (void) delta_min;

    if (delta_max != 1.0)
      max_eigv_squared *= delta_max;

    /* delta = sigma_min(a D_ker(0)), the leftmost real mode of the
    *unshifted* kernel (mass = 0, not the Kernel operator's -rho shift): this
    is the quantity the resonance condition is anchored to, not lambda_min
    of X. It is an ensemble constant (measured stable to ~1% across
    configurations), so if the caller already knows it, 'kernel_delta > 0'
    supplies it directly and skips the second Lanczos run entirely.
    Otherwise it is measured with a second, looser Lanczos run since only
    ~2 decimal places are needed, with an iteration floor since Lanczos
    approaches lambda_min from above and can plateau before a loose
    tolerance is met. */
    qpb_double delta;
    if(kernel_delta > 0.0)
    {
      delta = kernel_delta;
      print(" delta (leftmost real kernel mode, supplied) = %.6f\n", delta);
    }
    else
    {
      qpb_double delta_min_eigv_squared, delta_max_eigv_squared;
      int delta_Lanczos_iters = qpb_extreme_eigenvalues_of_X_squared(
                        &delta_min_eigv_squared, &delta_max_eigv_squared, 0.0,
                        1e-4, Lanczos_max_iters, 100);
      print(" Total number of Lanczos algorithm iterations (delta) = %d\n", \
                                                            delta_Lanczos_iters);
      delta = sqrt(delta_min_eigv_squared);
      print(" delta (leftmost real kernel mode)        = %.6f\n", delta);
      print(" sanity check: unshifted kernel max eigv^2 = %.6f "
            "(expect ~64-67 for rho=1, c_sw=0, Wilson kernel)\n", \
                                                          delta_max_eigv_squared);
    }

    qpb_double max_eigv = sqrt(max_eigv_squared);
    qpb_double lambda_min = sqrt(min_eigv_squared);   /* raw, pre-reduction */
    qpb_double theta_window = 0.5;   /* window/margin gate; see select_delta_min -
                                      unrelated to theta_gate below */
    qpb_double selected_delta_min = select_delta_min(min_eigv_squared, \
                      max_eigv, rho, mass, delta, Zol_iters, theta_window);
    if(selected_delta_min <= 0.0)
    {
      error(" !\n");
      error(" select_delta_min: no value of alpha (with delta_min <= 0.5) "
            "keeps the margin sigma_min >= theta_window*am across the "
            "window [rho - 1.1*delta, rho - 0.9*delta].\n");
      error(" This Zolotarev order (n = %d) cannot safely represent am = "
            "%g at rho = %g with delta = %g; raise n.\n", Zol_iters, mass, \
                                                                    rho, delta);
      error(" !\n");
      exit(QPB_PARAMETERS_ERROR);
    }

    min_eigv_squared *= selected_delta_min;

    print(" Min eigenvalue squared = %.16f\n", min_eigv_squared);
    print(" Max eigenvalue squared = %.16f\n", max_eigv_squared);
    print(" selected delta_min                       = %.6e\n", \
                                                            selected_delta_min);

    /* And then their square root value is stored inside the 'min_eigv' and
    'max_eigv' attributes of the 'ov_params' struct. */

    ov_params.min_eigv = sqrt(min_eigv_squared);
    ov_params.max_eigv = max_eigv;

    print(" resulting alpha = sqrt(delta_min*lam2min) = %.6e\n", \
                                                            ov_params.min_eigv);

    /* ----------------------- expansion coefficients ----------------------- */

    Zolotarev_order = Zol_iters;
    MS_solver_precision = ms_epsilon;
    MS_maximum_solver_iterations = ms_max_iter;

    /* Calculate the numerical terms of the partial fraction expansion */
    shifts = qpb_alloc(sizeof(qpb_double)*Zolotarev_order);
    numerators = qpb_alloc(sizeof(qpb_double)*Zolotarev_order);

    /* Allocate arrays */
    qpb_double *c = qpb_alloc(sizeof(qpb_double)*2*Zolotarev_order);
    qpb_double *b = qpb_alloc(sizeof(qpb_double)*Zolotarev_order);
    qpb_double normalization_constant;

    build_Zolotarev_coeffs(ov_params.min_eigv, ov_params.max_eigv, \
                Zolotarev_order, c, b, &normalization_constant);

    /* --------------------------- diagnostics --------------------------- */

    qpb_double C = rho + 0.5*mass, R = rho - 0.5*mass;
    qpb_double Z_at_rho_minus_delta = Zolotarev_at(rho - delta, \
                ov_params.min_eigv, Zolotarev_order, c, normalization_constant);
    qpb_double sigma_min = C - R*Z_at_rho_minus_delta;
    qpb_double Z_at_max = Zolotarev_at(ov_params.max_eigv, ov_params.min_eigv, \
                Zolotarev_order, c, normalization_constant);
    qpb_double kappa_CGNR = pow((C + R*Z_at_max)/sigma_min, 2);

    print(" Z_n(rho - delta)                         = %.6f\n", \
                                                          Z_at_rho_minus_delta);
    print(" predicted sigma_min(D_ov)                = %.6e\n", sigma_min);
    if(mass != 0.0)
      print(" sigma_min / am                           = %.6f\n", \
                                                                sigma_min/mass);
    print(" predicted kappa_CGNR                     = %.6e\n", kappa_CGNR);

    print_overshoot_regions(lambda_min, rho - delta, ov_params.min_eigv, C, R, \
                Zolotarev_order, c, normalization_constant);

    /* d_Z quality gate: the window/margin test above guards invertibility
    and index correctness near rho - delta, but says nothing about chiral
    symmetry breaking, which is governed by the equioscillation amplitude
    d_Z over the whole interval. theta_gate is a decided constant (0.75),
    not a per-run default, and a FAIL is a report, not an abort: low-order
    runs are legitimate convergence-study points, just not physics
    measurements (see ZOLOTAREV_PATCH_3.md). */
    qpb_double d_Z = equioscillation_amplitude(ov_params.min_eigv, \
                ov_params.max_eigv, Zolotarev_order);
    qpb_double theta_gate = 0.75;
    qpb_double d_Z_gate = theta_gate*mass /
                      (rho - mass/2.0 + theta_gate*mass/2.0);
    int chiral_quality_pass = (d_Z <= d_Z_gate);

    print(" d_Z (equioscillation amplitude)          = %.6e\n", d_Z);
    print(" d_Z gate  theta_gate*am/(rho-am/2+...)    = %.6e\n", d_Z_gate);
    print(" chiral quality                           = %s\n", \
                                    chiral_quality_pass ? "PASS" : "FAIL");
    if(!chiral_quality_pass)
      print(" WARNING: d_Z exceeds the gate. Index is protected by the "
            "window\n"
            "          guard, but chiral symmetry breaking is O(d_Z), "
            "which may\n"
            "          be comparable to or larger than am. Treat this run "
            "as a\n"
            "          convergence-study point, not a physics "
            "measurement.\n");

    constant_term = normalization_constant / ov_params.min_eigv;

    for(int i=0; i<Zolotarev_order; i++)
    {
      shifts[i] = c[2*i] * ov_params.min_eigv * ov_params.min_eigv;
      numerators[i] = normalization_constant * b[i] * ov_params.min_eigv;
    }

    free(c);
    free(b);

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
