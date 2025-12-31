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


#define OVERLAP_NUMB_TEMP_VECS 8
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
      for (int k = 1; k <= n-1; k++)
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

/* ------------------------ MATRIX-VECTOR FUNCTIONS ------------------------ */

void
qpb_overlap_Zolotarev_init(void * gauge, qpb_clover_term clover, \
          int Zol_iters, qpb_double rho, \
          qpb_double c_sw, qpb_double mass, qpb_double scaling_factor, \
          qpb_double ms_epsilon, int ms_max_iter, \
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

    constant_term = c[2*Zolotarev_order-1] * ov_params.min_eigv * ov_params.min_eigv;

    for(int i=0; i<Zolotarev_order; i++)
    {
      shifts[i] = c[2*i] * ov_params.min_eigv * ov_params.min_eigv;
      numerators[i] = normalization_constant * b[i] / ov_params.min_eigv;
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


INLINE void
X_squared_op(qpb_spinor_field y, qpb_spinor_field x)
{
  /* Implements γ5( (a*D - ρ)( γ5( (a*D - ρ)(x) ) ) ) */

  qpb_spinor_field z = ov_temp_vecs[0];

  void *dslash_args[4];

  dslash_args[0] = ov_params.gauge_ptr;
  dslash_args[1] = &ov_params.m_bare;
  dslash_args[2] = &ov_params.clover;
  dslash_args[3] = &ov_params.c_sw;

  ov_params.g5_dslash_op(z, x, dslash_args);
  ov_params.g5_dslash_op(y, z, dslash_args);

  return;
}


void
qpb_gamma5_sign_function_of_X_Zolotarev(qpb_spinor_field y, qpb_spinor_field x)
{
  /* Implements: γ5(sign(X(x))) = γ5(X(c_0 + Sum_{i=1}^{n} c_i/(X^2+σ_i) )),
      with X(x) = γ5(D(x) - ρ*x) . */

  qpb_spinor_field sum = ov_temp_vecs[1];

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

  // Initialize sum to zero
  qpb_spinor_field_set_zero(sum);  

  // And then add the rest of the partial fraction terms
  for(int sigma=0; sigma<Zolotarev_order; sigma++)
    qpb_spinor_axpy(sum, (qpb_complex) {numerators[sigma], 0.}, yMS[sigma], sum);

  // Act with X^2 + c[2n] term
  X_squared_op(y, sum);
  qpb_spinor_axpy(sum, (qpb_complex) {constant_term, 0.}, sum, y);

  // Act with gamma5 X term
  D_op(y, sum);

  return;
}


void
qpb_overlap_Zolotarev(qpb_spinor_field y, qpb_spinor_field x)
{
  /* Implements:
        Dov,m(x) = (rho+overlap_mass/2)*x+ (rho-overlap_mass/2)*g5(sign(X))
  */
  
  qpb_spinor_field z = ov_temp_vecs[2];

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
  qpb_spinor_field p = ov_temp_vecs[3];
  qpb_spinor_field r = ov_temp_vecs[4];
  qpb_spinor_field y = ov_temp_vecs[5];
  qpb_spinor_field w = ov_temp_vecs[6];
  qpb_spinor_field bprime = ov_temp_vecs[7];

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
