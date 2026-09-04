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
#include <qpb.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>

#define OVERLAP_NUMB_TEMP_VECS 7
#define MSCG_NUMB_TEMP_VECS 20


static qpb_spinor_field ov_temp_vecs[OVERLAP_NUMB_TEMP_VECS];
static qpb_spinor_field mscg_temp_vecs[MSCG_NUMB_TEMP_VECS];

static qpb_overlap_params ov_params;

static int KL_diagonal_order;
static qpb_double MS_solver_precision;
static int MS_maximum_solver_iterations;

static qpb_double *numerators;
static qpb_double *shifts;
static qpb_double constant_term;


/* --------------------- EXTREME EIGENVALUES FUNCTIONS --------------------- */

static void
tridiag_eigenv(double *eig, double *a, double *b, int n)
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


static int
qpb_extreme_eigenvalues_of_X_squared(qpb_double *min_eigv, \
  qpb_double *max_eigv, qpb_double bare_mass, qpb_double Lanczos_epsilon, \
  int max_iters, int min_iters)
{
  /* It calculates the extreme eigenvalues of the eigenvalue spectrum
  of H^2, H ≡ γ5*Kernel(x), with: Kernel(x) = (a*D + bare_mass)(x), using the
  Lanczos algorithm.

  'bare_mass' selects which kernel is probed: pass ov_params.m_bare = -rho
  for the shifted kernel operator X = γ5*(a*D - rho), whose spectrum is the
  one the sign function approximation has to cover.

  'min_iters' is a floor on the number of iterations before the convergence
  check is allowed to stop the loop: Lanczos approaches lambda_min from above
  and can sit on a plateau, so a loose tolerance combined with a small
  iteration count can stop early on a not-yet-converged estimate (pass 0 for
  no floor). */

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
  int n_ritz = 0;  /* size of the last tridiagonal system actually solved */
  for(i=1; i<max_iters; i++)
  {
    qpb_lanczos(a, b, solver_arg_links, clover_term, kappa, c_sw, -1);
    tridiag_eigenv(eig, a, b, i+1);
    n_ritz = i+1;

    lambda = eig[i] / eig[0];
    dlambda = fabs(lambda - lambda0) / fabs(lambda + lambda0);
    if (i%100==0)
      print("\titer = %4d, CN = %e/%e = %e (change = %e, target = %e)\n", i+1,\
                      eig[i], eig[0], eig[i]/eig[0], dlambda, Lanczos_epsilon);
    if((i+1 >= min_iters) && dlambda < Lanczos_epsilon*0.5)
      break;
    lambda0 = lambda;
  }

  /* The last tridiag_eigenv() call filled eig[0..n_ritz-1], so the largest
  Ritz value is eig[n_ritz-1]: indexing it as eig[i-1] is only correct when
  the loop ran to exhaustion, and returns the *second* largest Ritz value on
  the (usual) converged-and-broke path - inconsistent with the convergence
  check just above, which already treats eig[i] as lambda_max. */
  *min_eigv = (qpb_double) eig[0];
  *max_eigv = (qpb_double) eig[n_ritz-1];

  free(a);
  free(b);
  free(eig);

  qpb_lanczos_finalize();

  return i;
}

/* ------------------------ MATRIX-VECTOR FUNCTIONS ------------------------ */

void
qpb_overlap_Neuberger_init(void * gauge, qpb_clover_term clover, \
          enum qpb_kl_classes kl_class, int kl_iters, qpb_double rho, \
          qpb_double c_sw, qpb_double mass, qpb_double scaling_factor, \
          qpb_double ms_epsilon, int ms_max_iter, int optimal_scaling, \
          qpb_double Lanczos_epsilon, int Lanczos_max_iters, \
          qpb_double delta_min, qpb_double delta_max)
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

    /* ------------------------ the scaling parameter ------------------------ */

    /* mu is the scaling parameter that rescales the argument of the sign
    function approximation as X -> X/sqrt(mu). It is either taken verbatim
    from the caller ('scaling_factor', i.e. the "Scaling factor" input file
    entry) or, when 'optimal_scaling' is nonzero, computed from the spectrum
    of X as the so-called optimal scaling

      mu = alpha*beta,  alpha = min|X|, beta = max|X|,  |X| = sqrt(X^2)

    The Neuberger approximation built below is

      h_n(x) = x*(1/n)*sum_i sec^2(theta_i)/(x^2 + tan^2(theta_i))
             = [(1+x)^(2n) - (1-x)^(2n)] / [(1+x)^(2n) + (1-x)^(2n)]

    which is invariant under x -> 1/x, since (x-1)^(2n) = (1-x)^(2n). The
    approximation error at the two ends of the rescaled spectrum,
    alpha/sqrt(mu) and beta/sqrt(mu), is therefore equal exactly when those
    two are reciprocals of each other, i.e. when mu = alpha*beta: that choice
    equalizes the error at the spectral endpoints and so minimizes the worst
    error over the whole spectrum. */

    qpb_double mu = scaling_factor;

    if(optimal_scaling)
    {
      /* The extrema of the eigenvalue spectrum of X^2, X = g5*(a*D - rho),
      are calculated and stored in 'min_eigv_squared' and 'max_eigv_squared'
      correspondingly. */
      qpb_double min_eigv_squared, max_eigv_squared;
      int Lanczos_iters = qpb_extreme_eigenvalues_of_X_squared( \
                      &min_eigv_squared, &max_eigv_squared, \
                      ov_params.m_bare,  // shifted kernel: bare mass = -rho
                      Lanczos_epsilon, Lanczos_max_iters, 100);
      print(" Total number of Lanczos algorithm iterations = %d\n", \
                                                                Lanczos_iters);
      print(" Min eigenvalue squared (raw) = %.16f\n", min_eigv_squared);
      print(" Max eigenvalue squared (raw) = %.16f\n", max_eigv_squared);

      /* Safeguards: Lanczos returns interior estimates - it approaches
      lambda_min from above and lambda_max from below - so the measured
      interval is widened before it is used, by shrinking the lower end
      (delta_min < 1) and stretching the upper end (delta_max > 1). */
      if (delta_min != 1.0)
        min_eigv_squared *= delta_min;
      if (delta_max != 1.0)
        max_eigv_squared *= delta_max;

      if(min_eigv_squared <= 0.0 || max_eigv_squared <= 0.0)
      {
        error(" !\n");
        error(" qpb_overlap_Neuberger_init: non-positive eigenvalue estimate "
              "of X^2 (min = %g, max = %g) after applying delta_min = %g, "
              "delta_max = %g.\n", min_eigv_squared, max_eigv_squared, \
                                                        delta_min, delta_max);
        error(" !\n");
        exit(QPB_PARAMETERS_ERROR);
      }

      /* And then their square roots, the extrema of |X| = sqrt(X^2), are
      stored inside the 'min_eigv' and 'max_eigv' attributes of the
      'ov_params' struct. */
      ov_params.min_eigv = sqrt(min_eigv_squared);
      ov_params.max_eigv = sqrt(max_eigv_squared);

      qpb_double alpha = ov_params.min_eigv;
      qpb_double beta = ov_params.max_eigv;

      mu = alpha*beta;

      print(" alpha = min|X|               = %.16f\n", alpha);
      print(" beta  = max|X|               = %.16f\n", beta);
      print(" condition number beta/alpha  = %.16f\n", beta/alpha);
      print(" optimal scaling mu = alpha*beta = %.16f\n", mu);
      if(scaling_factor != mu)
        print(" NOTE: optimal scaling is on, so the supplied scaling factor "
              "(%g) is ignored.\n", scaling_factor);
    }

    /* ----------------------- expansion coefficients ----------------------- */

    KL_diagonal_order = kl_iters;
    MS_solver_precision = ms_epsilon;
    MS_maximum_solver_iterations = ms_max_iter;

    /* Calculate the numerical terms of the partial fraction expansion.
    'constant_term' is scratch storage for 1/n here; it is NOT an additive
    term of the expansion, and is not read again after this loop. */
    shifts = qpb_alloc(sizeof(qpb_double)*KL_diagonal_order);
    numerators = qpb_alloc(sizeof(qpb_double)*KL_diagonal_order);

    constant_term = 1.0/((qpb_double) (KL_diagonal_order));

    for(int i=0; i<KL_diagonal_order; i++)
    {
      qpb_double trig_arg = 0.5*M_PI*(i+0.5)*constant_term;
      shifts[i] = pow(tan(trig_arg), 2);
      numerators[i] = constant_term/powl(cos(trig_arg), 2);
      // print("numerator[%d] = %.25f, shift[%d] = %.25f\n", i, numerators[i], \
                                                              i, shifts[i]);
    }

    /* Modify the numerical constants of the partial fraction expansion using
    the scaling parameter mu, which rescales the argument of the
    approximation as X -> X/sqrt(mu):

      h_n(x; mu) = h_n(x/sqrt(mu))
                 = x*( c_0/sqrt(mu) + sum_i b_i*sqrt(mu)/(x^2 + mu*s_i) )

    i.e. b_i -> b_i*sqrt(mu), s_i -> s_i*mu. Note there is no additive
    constant to rescale: unlike the diagonal KL expansion, this one is a
    proper rational function with no c_0 term (see
    qpb_gamma5_sign_function_of_X_Neuberger below), and 'constant_term' here
    is only scratch storage holding 1/n while the coefficients are built. */
    if (mu != 1.0)
    {
      for(int i=0; i<KL_diagonal_order; i++)
      {
        numerators[i] *= sqrt(mu);
        shifts[i] *= mu;
      }
    }

    qpb_mscongrad_init(KL_diagonal_order);

  }
	
  return;
}


void
qpb_overlap_Neuberger_finalize()
{
  qpb_comm_halo_spinor_field_finalize();
  for(int i=0; i<OVERLAP_NUMB_TEMP_VECS; i++)
    qpb_spinor_field_finalize(ov_temp_vecs[i]);
  
  for(int i=0; i<MSCG_NUMB_TEMP_VECS; i++)
    qpb_spinor_field_finalize(mscg_temp_vecs[i]);

  if(which_dslash_op == QPB_DSLASH_STANDARD)
    qpb_gauge_field_finalize(*(qpb_gauge_field *)ov_params.gauge_ptr);
  
  ov_params.initialized = 0;
  
  qpb_mscongrad_finalize(KL_diagonal_order);

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
qpb_gamma5_sign_function_of_X_Neuberger(qpb_spinor_field y, qpb_spinor_field x)
{
  /* Implements: γ5(sign(X(x))) = γ5(X(Sum_{i=1}^{n} c_i/(X^2+σ_i))),
      with X(x) = γ5(D(x) - ρ*x) .

  Note there is no additive c_0 term here, unlike the diagonal KL expansion:
  with c_i = (1/n)*sec^2(θ_i) and σ_i = tan^2(θ_i), θ_i = (π/2n)*(i+1/2),
  this sum is exactly the Neuberger rational function
  [(1+x)^(2n) - (1-x)^(2n)] / [(1+x)^(2n) + (1-x)^(2n)], which already
  satisfies h_n(1) = 1. Adding a c_0 would give h_n(1) = (n+1)/n instead. */

  qpb_spinor_field sum = ov_temp_vecs[0];

  qpb_spinor_field yMS[KL_diagonal_order];
  for(int i=0; i<KL_diagonal_order; i++)
  {
    yMS[i] = mscg_temp_vecs[i];
    qpb_spinor_field_set_zero(yMS[i]);
  }

  qpb_double kernel_mass = ov_params.m_bare; // Kernel operator bare mass
  qpb_double kernel_kappa = 1./(2*kernel_mass+8.);

  qpb_mscongrad(yMS, x, ov_params.gauge_ptr, ov_params.clover, kernel_kappa, \
    KL_diagonal_order, shifts, ov_params.c_sw, MS_solver_precision, \
    MS_maximum_solver_iterations);

  // Add the partial fraction terms
  qpb_spinor_field_set_zero(sum);
  for(int i=0; i<KL_diagonal_order; i++)
    qpb_spinor_axpy(sum, (qpb_complex) {numerators[i], 0.}, yMS[i], sum);

  D_op(y, sum);

  return;
}


void
qpb_overlap_Neuberger(qpb_spinor_field y, qpb_spinor_field x)
{
  /* Implements:
        Dov,m(x) = (rho+overlap_mass/2)*x+ (rho-overlap_mass/2)*g5(sign(X))
  */
  
  qpb_spinor_field z = ov_temp_vecs[1];

  qpb_double overlap_mass = ov_params.mass; // Overlap operator Dov,m mass
  qpb_double rho = ov_params.rho;

  qpb_complex a = {rho + 0.5*overlap_mass, 0.};
  qpb_complex b = {rho - 0.5*overlap_mass, 0.};

  qpb_gamma5_sign_function_of_X_Neuberger(z, x);

  qpb_spinor_axpby(y, a, x, b, z);

  return;
}


void
qpb_gamma5_overlap_Neuberger(qpb_spinor_field y, qpb_spinor_field x)
{
  qpb_overlap_Neuberger(y, x);
  qpb_spinor_gamma5(y, y);

  return;
}


int
qpb_congrad_overlap_Neuberger(qpb_spinor_field x, qpb_spinor_field b, \
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
  qpb_gamma5_overlap_Neuberger(bprime, w);

  qpb_spinor_xdotx(&b_norm, bprime);

  qpb_spinor_field_set_zero(x);

  /* r0 = bprime - A(x) */
  // qpb_gamma5_overlap_Neuberger(w, x);
  // qpb_gamma5_overlap_Neuberger(p, w);
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
    qpb_gamma5_overlap_Neuberger(w, p);
    qpb_gamma5_overlap_Neuberger(y, w);

    /* omega = dot(p, A(p)) */
    qpb_spinor_xdoty(&omega, p, y);

    /* alpha = dot(r, r)/omega */
    alpha = CDEV(gamma, omega);

    /* x <- x + alpha*p */
    qpb_spinor_axpy(x, alpha, p, x);

    if(iters % n_reeval == 0) 
    {
      qpb_gamma5_overlap_Neuberger(w, x);
      qpb_gamma5_overlap_Neuberger(y, w);
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

  qpb_gamma5_overlap_Neuberger(w, x);
  qpb_gamma5_overlap_Neuberger(y, w);
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
