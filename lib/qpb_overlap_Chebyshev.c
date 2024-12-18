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
#include <math.h>
#include <qpb.h>

#include <qpb_overlap_Chebyshev_defs.h>
#include <qpb_overlap_Chebyshev.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>

#define OVERLAP_NUMB_TEMP_VECS 21

static qpb_spinor_field ov_temp_vecs[OVERLAP_NUMB_TEMP_VECS];

static int number_of_Chebyshev_terms;

static qpb_double *discrete_r_values;
static qpb_double *discrete_xk;
static qpb_double *Chebyshev_TnM2_for_k_array;
static qpb_double *Chebyshev_TnM1_for_k_array;

static qpb_complex *Chebyshev_coeff;

static qpb_spinor_field Chebyshev_Tn_vec;
static qpb_spinor_field Chebyshev_Tnm1_vec;

static qpb_Chebyshev_overlap_params ov_params;


/* -------------- SCALAR CHEBYSHEV POLYNOMIAL TERMS FUNCTIONS -------------- */

INLINE qpb_double
scalar_r_func(qpb_double x)
{
  /* It implements the r(x) function, defined as:

        r(x) = [ (beta^2 + alpha^2)/2 + x*(beta^2 - alpha^2)/2 ]^(-1/2),
  
     with alpha and beta the minimum and maximum (real) eigenvalues of the
     eigenvalues spectrum of H^2.
  */

  qpb_double alpha = ov_params.min_eigv;
  qpb_double beta = ov_params.max_eigv;

  qpb_double denominator = sqrt(0.5*(beta*beta + alpha*alpha) \
                                            + 0.5*x*(beta*beta - alpha*alpha));

  return 1./denominator;
}


INLINE qpb_double
scalar_nth_Chebyshev_polynomial_term_for_xk(int k, int n)
{
  /* Implementation of the recursive relation for the calculation of the nth
     degree (scalar) Chebyshev polynomial term of the 1st kind using the two
     terms of immediately previous degrees:
        
                      T_n(x) = 2*x*T_{n-1}(x) - T_{n-2}(x),
      
     with T_{n=0}(x) = 1 and T_{n=1}(x) = x, specifically for
     x = xk ≡ cos[ π/N*(k + 1/2) ], for k=0,..., N-1.

    NOTE 1: This function has access to the values of the Chebyshev polynomial
    terms of n-1 and n-2 degrees (for x=xk), T_{n-1}(x=xk) and T_{n-2}(x=xk),
    via the stored values inside the 'Chebyshev_TnM1_for_k_array' and 
    'Chebyshev_TnM2_for_k_array' arrays correspondingly for the corresponding k.

    Note 2: At the end, it stores the calculated nth degree Chebyshev polynomial
    term for x=xk inside the 'Chebyshev_TnM1_for_k_array' array (at index k),
    while previously the 'Chebyshev_TnM1_for_k_array'[k] value itself has been
    saved inside 'Chebyshev_TnM2_for_k_array'[k] for future use.

    Note 3: Furthermore, the function is called inside the loop for the
    calculation of the nth expansion coefficient c_n of the weighted expansion
    sum of the inverse square root of H^2. The expansion coefficients are
    calculated serially, from n=0 to n=N-1 with step 1. Therefore no need to
    check if the Chebyshev polynomial terms of the two previous degrees (for
    x=xk) have been already calculated, they were indeed by construction.
  */

  qpb_double Tn, TnM1, TnM2, xk;

  TnM2 = Chebyshev_TnM2_for_k_array[k];
  TnM1 = Chebyshev_TnM1_for_k_array[k];
  xk = discrete_xk[k];

  Tn = 2*xk*TnM1 - TnM2;

  Chebyshev_TnM2_for_k_array[k] = Chebyshev_TnM1_for_k_array[k];
  Chebyshev_TnM1_for_k_array[k] = Tn;

  return Tn;
}


qpb_complex
expansion_Chebyshev_coeff_n(int n)
{
  /* It calculates the expansion coefficient c_n for the nth term of the
     truncated series expansion of *matrix* Chebyshev polynomials approximating
     the inverse square root of H^2. The expression for the coefficient c_n
     is itself a truncated series of *scalar* Chebyshev polynomials:
 
              c_n = 1/(p_n N) sum_{k=0}^{N-1} r(x_k)*T_n(x_k),
 
     for x_k = cos[ π/N*(k + 1/2) ], and p_0 = 1, p_n = 1/2 for n>0.
 
     Note: This function assumes that x_k and r(x_k) values, for k = 1, ..., N,
     have been already calculated and stored inside the 'discrete_xk' and 
     'discrete_r_values' arrays correspondingly.
  */
  
  qpb_double sum_prefactor = 2./(qpb_double) number_of_Chebyshev_terms;

  qpb_double sum = 0.;
  if (n == 0)
  {
    sum_prefactor *= 0.5; // Necessary correction for the n=0 prefactor
    for (int k=0; k<number_of_Chebyshev_terms; ++k)
      sum += discrete_r_values[k]; // Simplified expression since T_{n=0}(x_k)=1
  }
  else if (n == 1)
    for (int k=0; k<number_of_Chebyshev_terms; ++k)
        sum += discrete_r_values[k]*discrete_xk[k]; // T_{n=1}(x_k)=x_k
  else
    for (int k=0; k<number_of_Chebyshev_terms; ++k)
      sum += discrete_r_values[k]\
                          *scalar_nth_Chebyshev_polynomial_term_for_xk(k, n);

  qpb_complex complex_Chebyshev_coefficient = {sum_prefactor*sum, 0};

  return complex_Chebyshev_coefficient;
}


/* --------------------- EXTREME EIGENVALUES FUNCTIONS --------------------- */

INLINE void
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
qpb_overlap_Chebyshev_init(void *gauge, qpb_clover_term clover, qpb_double rho,\
    qpb_double c_sw, qpb_double mass, qpb_double Lanczos_epsilon, \
      int Lanczos_max_iters, int N_Cheb, \
        qpb_double delta_max, qpb_double delta_min)
{
  if(ov_params.initialized != QPB_OVERLAP_INITIALIZED)
  {
    qpb_comm_halo_spinor_field_init();

    /* 'qpb_spinor_field' variables for storing already calculated (vector)
    Chebyshev polynomial terms of degree n and n-1 correspondingly, to be used
    in the calculation of the next-degree (n+1) term. They are initialized
    appropriately when function 'qpb_gamma5_sign_function_of_H' is called. */
    Chebyshev_Tn_vec = qpb_spinor_field_init();
    Chebyshev_Tnm1_vec = qpb_spinor_field_init();

    for(int i=0; i<OVERLAP_NUMB_TEMP_VECS; i++)
    {
      ov_temp_vecs[i] = qpb_spinor_field_init();
      qpb_spinor_field_set_zero(ov_temp_vecs[i]);
    }

    qpb_gauge_field gauge_bc;
    if(which_dslash_op == QPB_DSLASH_STANDARD)
    {
      gauge_bc = qpb_gauge_field_init();
      qpb_timebc_set_gauge_field(gauge_bc, *(qpb_gauge_field *)gauge, \
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
    ov_params.mass = mass; // Overlap operator mass
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


    /* --------------------- extreme eigenvalues of H^2 --------------------- */

    qpb_double min_eigv_squared;
    qpb_double max_eigv_squared;

    /* First the the extrema of the eigenvalues spectrum of H^2,
    H = g5*(D - rho), are calculated and are stored inside the
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

    number_of_Chebyshev_terms = N_Cheb;

    // Allocate arrays for caching calculated results
    discrete_r_values = qpb_alloc(sizeof(qpb_double)*number_of_Chebyshev_terms);
    discrete_xk = qpb_alloc(sizeof(qpb_double)*number_of_Chebyshev_terms);
    Chebyshev_TnM2_for_k_array = qpb_alloc(sizeof(qpb_double)\
                                                    *number_of_Chebyshev_terms);
    Chebyshev_TnM1_for_k_array = qpb_alloc(sizeof(qpb_double)\
                                                    *number_of_Chebyshev_terms);

    // Initialize all the caching arrays
    for(int k=0; k<number_of_Chebyshev_terms; k++)
    {
      qpb_double factor = M_PI/(qpb_double) number_of_Chebyshev_terms;
      qpb_double xk = cos(factor*(k+0.5));
      discrete_xk[k] = xk;
      discrete_r_values[k] = scalar_r_func(xk);
      Chebyshev_TnM2_for_k_array[k] = 1.;
      Chebyshev_TnM1_for_k_array[k] = discrete_xk[k];
    }

    // Allocate array for storing the expansion coefficients values
    Chebyshev_coeff = qpb_alloc(sizeof(qpb_complex)*number_of_Chebyshev_terms);
    /* Calculate the expansion coefficients of the weighted sum of the inverse
    square root of H^2. */
    for(int n=0; n<number_of_Chebyshev_terms; n++)
      Chebyshev_coeff[n] = expansion_Chebyshev_coeff_n(n);
    // // Print the expansion coefficients 
    // for(int n=0; n<number_of_Chebyshev_terms; n++)
    //   print("Chebyshev_coeff[%d]: %.12f\n", n, Chebyshev_coeff[n]);
  }

  return;
}


void
qpb_overlap_Chebyshev_finalize()
{
  qpb_comm_halo_spinor_field_finalize();

  for(int i=0; i<OVERLAP_NUMB_TEMP_VECS; i++)
    qpb_spinor_field_finalize(ov_temp_vecs[i]);

  qpb_spinor_field_finalize(Chebyshev_Tn_vec);
  qpb_spinor_field_finalize(Chebyshev_Tnm1_vec);

  free(Chebyshev_coeff);
  free(Chebyshev_TnM2_for_k_array);
  free(Chebyshev_TnM1_for_k_array);
  free(discrete_r_values);
  free(discrete_xk);
  
  if(which_dslash_op == QPB_DSLASH_STANDARD)
  {
    qpb_gauge_field_finalize(*(qpb_gauge_field *)ov_params.gauge_ptr);
  }

  qpb_lanczos_finalize();

  ov_params.initialized=0;

  return;
}


INLINE void
D_op(qpb_spinor_field y, qpb_spinor_field x)
{
  /* Kernel(x) = (a*D - ρ)(x) */

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
  /* H(x) ≡ γ5( (a*D - ρ)(x) ) */

  void *dslash_args[4];

  dslash_args[0] = ov_params.gauge_ptr;
  dslash_args[1] = &ov_params.m_bare;
  dslash_args[2] = &ov_params.clover;
  dslash_args[3] = &ov_params.c_sw;

  ov_params.g5_dslash_op(y, x, dslash_args);

  return;
}


INLINE void
Q_op(qpb_spinor_field y, qpb_spinor_field x)
{
  /* It implements:
     
        X(x) = ( 2*H(H(x)) - (beta^2+alpha^2)*x )/(beta^2-alpha^2)
             = 2/(beta^2-alpha^2)*H(H(x)) - (beta^2+alpha^2)/(beta^2-alpha^2)*x
             ≡ a*H(H(x)) - b*x,
     
     with alpha and beta the minimum and maximum (real) eigenvalues of the
     eigenvalues spectrum of H^2.
  */

  qpb_spinor_field z = ov_temp_vecs[0];

  qpb_double alpha = ov_params.min_eigv;
  qpb_double beta = ov_params.max_eigv;

  qpb_double prefactor = 1./(beta*beta - alpha*alpha);
  qpb_complex a = {2*prefactor, 0.0};
  qpb_complex b = {-prefactor*(beta*beta + alpha*alpha), 0.0};
  
  X_op(y, x);
  X_op(z, y);
  
  qpb_spinor_axpby(y, a, z, b, x);

  return;
}


INLINE void
qpb_Tnp1_Chebyshev_polynomial_term(qpb_spinor_field y)
{
  /* Implementation of the recursive relation for the calculation of the
     nth degree matrix Chebyshev polynomial term of the 1st kind, with the
     use of the two Chebyshev polynomials of immediately previous degrees:
 
                  T_n(x) = 2*X(T_{n-1}(x)) - T_{n-2}(x),
 
     with T_{n=0}(X(x)) = x, T_{n=1}(X(x)) = X(x).
 
     Note 1: At the end, the calculated T_n(x) vector is stored inside the 
     'Chebyshev_Tn_vec' 'qpb_spinor_field' variable, while previously the 
     T_{n-1}(x) vector (from 'Chebyshev_Tn_vec') is stored once again inside
     'Chebyshev_Tnm1_vec' this time, for future use.
 
     NOTE 2: This function is called inside the loop for the calculation of the
     weighted expansion sum of the inverse square root of H^2.
     This for loop repeats with increasing n from 0 to N-1 with step 1.
     Therefore no need to check if the Chebyshev polynomial terms of the two
     previous degrees have been already calculated, they were indeed by
     construction.
  */

  qpb_spinor_field z = ov_temp_vecs[1];

  qpb_complex two = {2., 0.};
  Q_op(y, Chebyshev_Tn_vec);
  qpb_spinor_ax(z, two, y);
  qpb_spinor_xmy(y, z, Chebyshev_Tnm1_vec);

  // Update Chebyshev_Tnm1_vec and Chebyshev_Tn_vec
  qpb_spinor_xeqy(Chebyshev_Tnm1_vec, Chebyshev_Tn_vec);
  qpb_spinor_xeqy(Chebyshev_Tn_vec, y);

  return;
}


void
qpb_gamma5_sign_function_of_H_for_N3(qpb_spinor_field y, qpb_spinor_field x)
{
  qpb_spinor_field sum = ov_temp_vecs[10];
  qpb_spinor_field z = ov_temp_vecs[11];

  qpb_double alpha = ov_params.min_eigv;
  qpb_double beta = ov_params.max_eigv;

  qpb_double prefactor = 1./(beta*beta - alpha*alpha);
  qpb_complex chi = {2*prefactor, 0.0};
  qpb_complex phi = {-prefactor*(beta*beta + alpha*alpha), 0.0};
  
  qpb_double a0 = Chebyshev_coeff[0].re - Chebyshev_coeff[2].re + Chebyshev_coeff[1].re*phi.re + 2*Chebyshev_coeff[2].re*phi.re*phi.re;
  qpb_double a1 = Chebyshev_coeff[1].re*chi.re + 4*Chebyshev_coeff[2].re*phi.re*chi.re;
  qpb_double a2 = 2*Chebyshev_coeff[2].re*chi.re*chi.re;

  print("a0=%.6f, a1=%.6f, a2=%.6f\n", a0, a1, a2);

  qpb_spinor_ax(sum, (qpb_complex) {a0, 0.0}, x);

  // X^2
  X_op(y, x);
  X_op(z, y);

  qpb_spinor_axpy(sum, (qpb_complex) {a1, 0.0}, z, sum);

  // X^4
  X_op(y, z);
  X_op(z, y);

  qpb_spinor_axpy(sum, (qpb_complex) {a2, 0.0}, z, sum);

  D_op(y, sum);

  return;
}


void
qpb_gamma5_sign_function_of_H_for_N4(qpb_spinor_field y, qpb_spinor_field x)
{
  qpb_spinor_field sum = ov_temp_vecs[10];
  qpb_spinor_field z = ov_temp_vecs[11];

  qpb_double alpha = ov_params.min_eigv;
  qpb_double beta = ov_params.max_eigv;

  qpb_double prefactor = 1./(beta*beta - alpha*alpha);
  qpb_complex chi = {2*prefactor, 0.0};
  qpb_complex phi = {-prefactor*(beta*beta + alpha*alpha), 0.0};
  
  qpb_double a0 = Chebyshev_coeff[0].re - Chebyshev_coeff[2].re + Chebyshev_coeff[1].re*phi.re - 3*Chebyshev_coeff[3].re*phi.re + 2*Chebyshev_coeff[2].re*phi.re*phi.re + 4*Chebyshev_coeff[3].re*phi.re*phi.re*phi.re;
  qpb_double a1 = Chebyshev_coeff[1].re*chi.re - 3*Chebyshev_coeff[3].re*chi.re + 4*Chebyshev_coeff[2].re*phi.re*chi.re + 12*Chebyshev_coeff[3].re*phi.re*phi.re*chi.re;
  qpb_double a2 = 2*Chebyshev_coeff[2].re*chi.re*chi.re + 12*Chebyshev_coeff[3].re*phi.re*chi.re*chi.re;
  qpb_double a3 = 4*Chebyshev_coeff[3].re*chi.re*chi.re*chi.re;

  print("a0=%.6f, a1=%.6f, a2=%.6f, a3=%.6f\n", a0, a1, a2, a3);

  qpb_spinor_ax(sum, (qpb_complex) {a0, 0.0}, x);

  // X^2
  X_op(y, x);
  X_op(z, y);

  qpb_spinor_axpy(sum, (qpb_complex) {a1, 0.0}, z, sum);

  // X^4
  X_op(y, z);
  X_op(z, y);

  qpb_spinor_axpy(sum, (qpb_complex) {a2, 0.0}, z, sum);

  // X^6
  X_op(y, z);
  X_op(z, y);

  qpb_spinor_axpy(sum, (qpb_complex) {a3, 0.0}, z, sum);

  D_op(y, sum);

  return;
}


void
qpb_Clenshaw_sum(qpb_spinor_field y, qpb_spinor_field x)
// qpb_gamma5_sign_function_of_H_Clenshaw(qpb_spinor_field y, qpb_spinor_field x)
{
  qpb_spinor_field z = ov_temp_vecs[17];

  qpb_spinor_field b_kp1 = ov_temp_vecs[18];
  qpb_spinor_field b_kp2 = ov_temp_vecs[19];

  qpb_spinor_field_set_zero(b_kp2);
  qpb_spinor_field_set_zero(b_kp1);

  qpb_complex two = {2., 0.};

  for(int n=number_of_Chebyshev_terms-1; n>0; n--)
  {
    // 1st term
    qpb_spinor_ax(z, Chebyshev_coeff[n], x);
    // 2nd term
    Q_op(y, b_kp1);
    qpb_spinor_axpy(z, two, y, z);
    // 3rd term
    qpb_spinor_xmy(y, z, b_kp2);
    
    qpb_spinor_xeqy(b_kp2, b_kp1);
    qpb_spinor_xeqy(b_kp1, y);
  }

  // Final sum
  qpb_spinor_ax(z, Chebyshev_coeff[0], x);
  Q_op(y, b_kp1);
  qpb_spinor_xpy(z, z, y);
  qpb_spinor_xmy(y, z, b_kp2);

  // γ5(H(x)) corresponds simply to D(x)
  // D_op(y, z);  
  
  return;
}


void
qpb_gamma5_sign_function_of_X(qpb_spinor_field y, qpb_spinor_field x)
{
  /* Implements: γ5( sign( H(x) ) ) = γ5( H( Sum_{n=0}^{N-1} c_n*T_n(x) ) ). */

  qpb_spinor_field sum = ov_temp_vecs[2];
  qpb_spinor_field Tn = ov_temp_vecs[3];

  // Initialize Chebyshev_Tnm1_vec with T_{n=0}
  qpb_spinor_xeqy(Chebyshev_Tnm1_vec, x);
  // Initialize Chebyshev_Tn_vec with T_{n=1}
  Q_op(y, x);
  qpb_spinor_xeqy(Chebyshev_Tn_vec, y);

  // Initialize sum with the term: c_{n=0}*T_{n=0}
  qpb_spinor_ax(sum, Chebyshev_coeff[0], Chebyshev_Tnm1_vec);
  // Initialize Tn with T_{n=1}
  qpb_spinor_xeqy(Tn, Chebyshev_Tn_vec);
  for (int n=1; n<number_of_Chebyshev_terms-1; ++n)
  {
    qpb_spinor_axpy(sum, Chebyshev_coeff[n], Tn, sum);
    // Calculate the next Chebyshev polynomial term Tnp1
    qpb_Tnp1_Chebyshev_polynomial_term(Tn);
  }
  // Add the last term to the sum
  qpb_spinor_axpy(sum, Chebyshev_coeff[number_of_Chebyshev_terms-1], Tn, sum);

  // γ5(H(x)) corresponds simply to D(x)
  D_op(y, sum);

  return;
}


void
qpb_overlap_Chebyshev(qpb_spinor_field y, qpb_spinor_field x)
{
  /* Implementing:
        Dov,m(x) = (rho+overlap_mass/2)*x+ (rho-overlap_mass/2)*g5(sign(X))
  */
  
  qpb_spinor_field z = ov_temp_vecs[4];

  qpb_double overlap_mass = ov_params.mass; // Overlap operator Dov,m mass
  qpb_double rho = ov_params.rho;

  qpb_complex a = {rho + 0.5*overlap_mass, 0.};
  qpb_complex b = {rho - 0.5*overlap_mass, 0.};

  // qpb_gamma5_sign_function_of_X(z, x);
  // qpb_gamma5_sign_function_of_H_Clenshaw(z, x);
  qpb_Clenshaw_sum(y, x);
  D_op(z, y);
  // qpb_gamma5_sign_function_of_H_for_N4(z, x);

  qpb_spinor_axpby(y, a, x, b, z);

  return;
}


void
qpb_overlap_Chebyshev_complex_conjugate(qpb_spinor_field y, qpb_spinor_field x)
{
  /* Implementing:
        Dov,m(x) = (rho+overlap_mass/2)*x+ (rho-overlap_mass/2)*g5(sign(X))
  */
  
  qpb_spinor_field z = ov_temp_vecs[4];

  qpb_double overlap_mass = ov_params.mass; // Overlap operator Dov,m mass
  qpb_double rho = ov_params.rho;

  qpb_complex a = {rho + 0.5*overlap_mass, 0.};
  qpb_complex b = {rho - 0.5*overlap_mass, 0.};

  // qpb_gamma5_sign_function_of_X(z, x);
  // qpb_gamma5_sign_function_of_H_Clenshaw(z, x);
  qpb_spinor_gamma5(z, x);
  X_op(y, z);
  qpb_Clenshaw_sum(z, y);
  // qpb_gamma5_sign_function_of_H_for_N4(z, x);

  qpb_spinor_axpby(y, a, x, b, z);

  return;
}


void
qpb_gamma5_overlap_Chebyshev(qpb_spinor_field y, qpb_spinor_field x)
{
  qpb_overlap_Chebyshev(y, x);
  
  qpb_spinor_gamma5(y, y);

  return;
}


int
qpb_congrad_overlap_Chebyshev(qpb_spinor_field x, qpb_spinor_field b,\
                                        qpb_double CG_epsilon, int CG_max_iter)
{
  qpb_spinor_field p = ov_temp_vecs[5];
  qpb_spinor_field r = ov_temp_vecs[6];
  qpb_spinor_field y = ov_temp_vecs[7];
  qpb_spinor_field w = ov_temp_vecs[8];
  qpb_spinor_field bprime = ov_temp_vecs[9];


  int n_reeval = 100;
  int n_echo = 10;
  int iters = 0;

  qpb_double res_norm, b_norm;
  qpb_complex_double alpha = {1, 0}, omega = {1, 0};
  qpb_complex_double beta, gamma;

  // bprime = Dov^+ b
  qpb_spinor_gamma5(w, b);
  qpb_gamma5_overlap_Chebyshev(bprime, w);

  qpb_spinor_xdotx(&b_norm, bprime);
  
  qpb_spinor_field_set_zero(x);

  /* r0 = bprime - A(x) */
  // qpb_gamma5_overlap_Chebyshev(w, x);
  // qpb_gamma5_overlap_Chebyshev(p, w);
  // qpb_spinor_xmy(r, bprime, p);

  /* Or r0=bprime for short since x0=0 */
  qpb_spinor_xeqy(r, bprime);

  qpb_spinor_xdotx(&gamma.re, r);
  gamma.im = 0;
  res_norm = gamma.re;
  /* p = r0 */
  qpb_spinor_xeqy(p, r);

  qpb_double t = qpb_stop_watch(0);
  for(iters=1; iters<CG_max_iter; iters++)
  {
    if(res_norm / b_norm <= CG_epsilon)
    {
      // print("CG stopped at: %.25e\n", res_norm / b_norm );
      break;
    }

    /* y = A(p) */
    qpb_gamma5_overlap_Chebyshev(w, p);
    qpb_gamma5_overlap_Chebyshev(y, w);

    /* omega = dot(p, A(p)) */
    qpb_spinor_xdoty(&omega, p, y);

    /* alpha = dot(r, r)/omega */
    alpha = CDEV(gamma, omega);

    /* x <- x + alpha*p */
    qpb_spinor_axpy(x, alpha, p, x);
    if(iters % n_reeval == 0) 
    {
      qpb_gamma5_overlap_Chebyshev(w, x);
      qpb_gamma5_overlap_Chebyshev(y, w);
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
	    print(" \t iters = %8d, res = %.25e\n", iters, res_norm / b_norm);

    beta.re = res_norm / gamma.re;
    beta.im = 0.;
    qpb_spinor_axpy(p, beta, p, r);
    gamma.re = res_norm;
    gamma.im = 0.;
  }

  t = qpb_stop_watch(t);
  
  // // Chebyshev inversion check
  // qpb_overlap_Chebyshev(y, x);
  // qpb_spinor_xmy(r, y, b);
  // qpb_double r_norm, b_original_norm;
  // qpb_spinor_xdotx(&r_norm, r);
  // qpb_spinor_xdotx(&b_original_norm, b);
  // print("Chebyshev inversion check: %.2e\n", r_norm/b_original_norm);

  qpb_gamma5_overlap_Chebyshev(w, x);
  qpb_gamma5_overlap_Chebyshev(y, w);
  qpb_spinor_xmy(r, bprime, y);
  qpb_spinor_xdotx(&res_norm, r);

  if(iters==CG_max_iter)
  {
    error(" !\n");
    error(" CG *did not* converge, after %d iterations\n", iters);
    error(" residual = %e, relative = %.25e, t = %g sec\n", res_norm, res_norm / b_norm, t);
    error(" !\n");
    return -1;
  }

  print(" \tAfter %d iters, CG converged, res = %e, relative = %e, "
        "t = %g sec\n",
         iters, res_norm, res_norm / b_norm, t);

  return iters;
}



int
qpb_congrad_overlap_Chebyshev_TEST(qpb_spinor_field x, qpb_spinor_field b,\
                                              qpb_double epsilon, int max_iter)
{  
  qpb_spinor_field p = ov_temp_vecs[5];
  qpb_spinor_field r = ov_temp_vecs[6];
  qpb_spinor_field y = ov_temp_vecs[7];
  qpb_spinor_field w = ov_temp_vecs[8];
  qpb_spinor_field bprime = ov_temp_vecs[9];


  int n_reeval = 100;
  int n_echo = 100;
  int iters = 0;

  qpb_double res_norm, b_norm;
  qpb_complex_double alpha = {1, 0}, omega = {1, 0};
  qpb_complex_double beta, gamma;

  // bprime = Dov^+ b
  qpb_overlap_Chebyshev_complex_conjugate(bprime, b);
  // qpb_spinor_gamma5(w, b);
  // qpb_gamma5_overlap_Chebyshev(bprime, w);

  qpb_spinor_xdotx(&b_norm, bprime);
  
  /* r0 = bprime - A(x) */
  qpb_overlap_Chebyshev(w, x);
  qpb_overlap_Chebyshev_complex_conjugate(p, w);
  qpb_spinor_xmy(r, bprime, p);

  qpb_spinor_xdotx(&gamma.re, r);
  gamma.im = 0;
  res_norm = gamma.re;
  /* p = r0 */
  qpb_spinor_xeqy(p, r);

  qpb_double t = qpb_stop_watch(0);
  for(iters=1; iters<max_iter; iters++)
  {
    if(res_norm / b_norm <= epsilon)
	    break;

    /* y = A(p) */
    qpb_overlap_Chebyshev(w, p);
    qpb_overlap_Chebyshev_complex_conjugate(y, w);

    /* omega = dot(p, A(p)) */
    qpb_spinor_xdoty(&omega, p, y);

    /* alpha = dot(r, r)/omega */
    alpha = CDEV(gamma, omega);

    /* x <- x + alpha*p */
    qpb_spinor_axpy(x, alpha, p, x);
    if(iters % n_reeval == 0) 
    {
      qpb_overlap_Chebyshev(w, x);
      qpb_overlap_Chebyshev_complex_conjugate(y, w);
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
  
  // Chebyshev inversion check
  qpb_overlap_Chebyshev(y, x);
  qpb_spinor_xmy(r, y, bprime);
  qpb_double r_norm, b_original_norm;
  qpb_spinor_xdotx(&r_norm, r);
  qpb_spinor_xdotx(&b_original_norm, b);
  print("Chebyshev inversion check: %.2e\n", r_norm/b_original_norm);

  qpb_overlap_Chebyshev(w, x);
  qpb_overlap_Chebyshev_complex_conjugate(y, w);
  // qpb_spinor_xmy(r, b, y);
  qpb_spinor_xmy(r, bprime, y);
  qpb_spinor_xdotx(&res_norm, r);

  if(iters==max_iter)
  {
    error(" !\n");
    error(" CG *did not* converge, after %d iterations\n", iters);
    error(" residual = %e, relative = %e, t = %g sec\n", res_norm,\
                                                        res_norm / b_norm, t);
    error(" !\n");
    return -1;
  }

  print(" \tAfter %d iters, CG converged, res = %e, relative = %e, "
        "t = %g sec\n",
         iters, res_norm, res_norm / b_original_norm, t);

  return iters;
}


int
qpb_bicgstab_overlap_Chebyshev(qpb_spinor_field x, qpb_spinor_field b,\
                                              qpb_double epsilon, int max_iter)
{
  qpb_spinor_field r0 = ov_temp_vecs[12];
  qpb_spinor_field r = ov_temp_vecs[13];
  qpb_spinor_field p = ov_temp_vecs[14];
  qpb_spinor_field u = ov_temp_vecs[15];
  qpb_spinor_field v = ov_temp_vecs[16];

  int iters = 0;
  const int n_reeval = 10000;
  int n_echo = 1;

  qpb_double res_norm, b_norm;
  qpb_complex_double alpha = {1, 0}, omega = {1, 0};
  qpb_complex_double beta, gamma, rho, zeta;
  // qpb_double mass = 1./(2.*kappa) - 4.;

  qpb_spinor_field_set_zero(x);

  qpb_spinor_xdotx(&b_norm, b);
  qpb_overlap_Chebyshev(r, x);
  qpb_spinor_xmy(r, b, r);
  qpb_spinor_xeqy(r0, r);
  qpb_spinor_xdotx(&gamma.re, r);
  gamma.im = 0;
  res_norm = gamma.re;
  rho = gamma;
  qpb_double t = qpb_stop_watch(0);
  for(iters=1; iters<max_iter; iters++)
    {
      if(res_norm / b_norm <= epsilon)
	break;
      
      qpb_spinor_xdoty(&gamma, r0, r);
      beta = CMUL(CDEV(gamma, rho), CDEV(alpha, omega));
      omega = CNEGATE(omega);
      qpb_spinor_axpy(p, omega, u, p);
      qpb_spinor_axpy(p, beta, p, r);
      qpb_overlap_Chebyshev(u, p);
      qpb_spinor_xdoty(&beta, r0, u);
      rho = gamma;
      alpha = CDEV(rho, beta);
      alpha = CNEGATE(alpha);
      qpb_spinor_axpy(r, alpha, u, r);
      qpb_overlap_Chebyshev(v, r);
      qpb_spinor_xdoty(&zeta, v, r);
      qpb_spinor_xdotx(&beta.re, v);
      beta.im = 0;
      omega = CDEV(zeta, beta);
      alpha = CNEGATE(alpha);
      qpb_spinor_axpy(x, omega, r, x);
      qpb_spinor_axpy(x, alpha, p, x);
      if(iters % n_reeval == 0)
	{
	  qpb_overlap_Chebyshev(r, x);
	  qpb_spinor_xmy(r, b, r);
	}
      else
	{
	  omega = CNEGATE(omega);
	  qpb_spinor_axpy(r, omega, v, r);
	  omega = CNEGATE(omega);
	}
      qpb_spinor_xdotx(&res_norm, r);
      
    if((iters % n_echo == 0))
	    print(" iters = %8d, res = %e\n", iters, res_norm / b_norm);
    }
  t = qpb_stop_watch(t);
  qpb_overlap_Chebyshev(r, x);
  qpb_spinor_xmy(r, b, r);
  qpb_spinor_xdotx(&res_norm, r);

  if(iters==max_iter)
    {
      error(" !\n");
      error(" BiCGStab *did not* converge, after %d iterations\n", iters);
      error(" residual = %e, relative = %e, t = %g sec\n", res_norm, res_norm / b_norm, t);
      error(" !\n");
      return -1;
    }

  print(" \tAfter %d iterations BiCGStab converged\n", iters);
  print(" residual = %e, relative = %e, t = %g sec\n", res_norm, res_norm / b_norm, t);
  
  return iters;
}

