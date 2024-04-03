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


#define OVERLAP_NUMB_TEMP_VECS 21
#define CG_NUMB_TEMP_VECS 10
#define MSCG_NUMB_TEMP_VECS 20


static qpb_spinor_field ov_temp_vecs[OVERLAP_NUMB_TEMP_VECS];
static qpb_spinor_field cg_temp_vecs[CG_NUMB_TEMP_VECS];
static qpb_spinor_field mscg_temp_vecs[MSCG_NUMB_TEMP_VECS];
static qpb_overlap_params ov_params;
static qpb_double scaling_factor;

void
qpb_overlap_kl_init(void * gauge, qpb_clover_term clover, qpb_double rho,\
                                qpb_double c_sw, qpb_double mass, qpb_double mu)
{
  if(ov_params.initialized != QPB_OVERLAP_INITIALIZED)
  {
    qpb_comm_halo_spinor_field_init();
    for(int i=0; i<OVERLAP_NUMB_TEMP_VECS; i++)
    {
      ov_temp_vecs[i] = qpb_spinor_field_init();
      qpb_spinor_field_set_zero(ov_temp_vecs[i]);
    }

    for(int i=0; i<CG_NUMB_TEMP_VECS; i++)
    {
      cg_temp_vecs[i] = qpb_spinor_field_init();
      qpb_spinor_field_set_zero(cg_temp_vecs[i]);
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

    scaling_factor = mu;
  }
	
  return;
}


void
qpb_overlap_kl_finalize()
{
  qpb_comm_halo_spinor_field_finalize();
  for(int i=0; i<OVERLAP_NUMB_TEMP_VECS; i++)
    qpb_spinor_field_finalize(ov_temp_vecs[i]);
  
  for(int i=0; i<CG_NUMB_TEMP_VECS; i++)
    qpb_spinor_field_finalize(cg_temp_vecs[i]);

  for(int i=0; i<MSCG_NUMB_TEMP_VECS; i++)
    qpb_spinor_field_finalize(mscg_temp_vecs[i]);

  if(which_dslash_op == QPB_DSLASH_STANDARD)
    qpb_gauge_field_finalize(*(qpb_gauge_field *)ov_params.gauge_ptr);
  
  ov_params.initialized = 0;
  
  return;
}


INLINE void
D_op(qpb_spinor_field y, qpb_spinor_field x)
{
  void *dslash_args[4];

  dslash_args[0] = ov_params.gauge_ptr;
  dslash_args[1] = &ov_params.m_bare;
  dslash_args[2] = &ov_params.clover;
  dslash_args[3] = &ov_params.c_sw;
  
  ov_params.dslash_op(y, x, dslash_args);
  
  return;
}


INLINE void
A_op(qpb_spinor_field z, qpb_spinor_field x)
{
  qpb_spinor_field y = ov_temp_vecs[0];
  
  void *dslash_args[4];
  
  dslash_args[0] = ov_params.gauge_ptr;
  dslash_args[1] = &ov_params.m_bare;
  dslash_args[2] = &ov_params.clover;
  dslash_args[3] = &ov_params.c_sw;

  ov_params.g5_dslash_op(z, x, dslash_args);
  ov_params.g5_dslash_op(y, z, dslash_args);
  
  qpb_complex three = {3.0,0.0};
  qpb_spinor_axpy(z, three, y, x);
  
  return;
}


INLINE void
XdaggerX_op(qpb_spinor_field z, qpb_spinor_field x)
{
  qpb_spinor_field y = ov_temp_vecs[1];
  
  void *dslash_args[4];
  
  dslash_args[0] = ov_params.gauge_ptr;
  dslash_args[1] = &ov_params.m_bare;
  dslash_args[2] = &ov_params.clover;
  dslash_args[3] = &ov_params.c_sw;

  ov_params.g5_dslash_op(y, x, dslash_args);
  ov_params.g5_dslash_op(z, y, dslash_args);
  
  return;
}


INLINE void
XdaggerX_plus_shift_op(qpb_spinor_field z, qpb_spinor_field x, qpb_double shift)
{
  qpb_spinor_field y = ov_temp_vecs[2];
  
  void *dslash_args[4];
  
  dslash_args[0] = ov_params.gauge_ptr;
  dslash_args[1] = &ov_params.m_bare;
  dslash_args[2] = &ov_params.clover;
  dslash_args[3] = &ov_params.c_sw;

  ov_params.g5_dslash_op(z, x, dslash_args);
  ov_params.g5_dslash_op(y, z, dslash_args);
  
  qpb_complex complex_shift = {shift, 0.0};
  qpb_spinor_axpy(z, complex_shift, x, y);

  return;
}


int
qpb_congrad_1p3A(qpb_spinor_field x, qpb_spinor_field b, qpb_double epsilon,\
                                                      int max_iter, int n_echo)
{  
  qpb_spinor_field p = cg_temp_vecs[0];
  qpb_spinor_field r = cg_temp_vecs[1];
  qpb_spinor_field y = cg_temp_vecs[2];

  int n_reeval = 100;
  int iters = 0;

  qpb_double res_norm, b_norm;
  qpb_complex_double alpha = {1, 0}, omega = {1, 0};
  qpb_complex_double beta, gamma;
  qpb_spinor_xdotx(&b_norm, b);
  
  /* r = (1+3A)x */
  A_op(p, x);
  /* r = b - (3*Ax + x) */
  qpb_spinor_xmy(r, b, p);
  qpb_spinor_xdotx(&gamma.re, r);
  gamma.im = 0;
  res_norm = gamma.re;
  qpb_spinor_xeqy(p, r);
  qpb_double t = qpb_stop_watch(0);
  for(iters=1; iters<max_iter; iters++)
  {
    if(res_norm / b_norm <= epsilon)
	    break;
    /* y = 3*Ap + p */
    A_op(y, p);

    /* omega = dot(p, (1+3*A)p) */
    qpb_spinor_xdoty(&omega, p, y);

    /* alpha = dot(r, r)/omega */
    alpha = CDEV(gamma, omega);

    /* x <- x + alpha*p */
    qpb_spinor_axpy(x, alpha, p, x);
    if(iters % n_reeval == 0) 
    {
      A_op(y, x);
      qpb_spinor_xmy(r, b, y);
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
  A_op(y, x);
  qpb_spinor_xmy(r, b, y);
  qpb_spinor_xdotx(&res_norm, r);

  if(iters==max_iter)
  {
    error(" !\n");
    error(" CG *did not* converge, after %d iterations\n", iters);
    error(" residual = %e, relative = %e, t = %g secs\n", res_norm,\
                                                        res_norm / b_norm, t);
    error(" !\n");
    return -1;
  }

  print(" \t After %d iters, CG converged, res = %e, relative = %e,\
                        t = %g secs\n", iters, res_norm, res_norm / b_norm, t);
  
  return iters;
}


void
qpb_gamma5_sign_function_partial_fractions(qpb_spinor_field y,\
    qpb_spinor_field x, qpb_complex constant_term, qpb_double* numerators,\
    qpb_double* shifts, int number_of_terms, qpb_double epsilon, int max_iter)
{
  qpb_spinor_field sum = ov_temp_vecs[3];
  qpb_spinor_field vector_diff = ov_temp_vecs[4];

  int n_echo = 100;

  qpb_mscongrad_init(number_of_terms);

  qpb_spinor_field yMS[number_of_terms];
  for(int sigma=0; sigma<number_of_terms; sigma++)
  {
    yMS[sigma] = mscg_temp_vecs[sigma];
    // It needs to re-initialized to 0 with every call of the function
    qpb_spinor_field_set_zero(yMS[sigma]);
  }

  qpb_double kernel_mass = ov_params.m_bare; // Kernel operator bare mass
  qpb_double kappa = 1./(2*kernel_mass+8.);

  qpb_mscongrad(yMS, x, ov_params.gauge_ptr, ov_params.clover, kappa,\
                  number_of_terms, shifts, ov_params.c_sw, epsilon, max_iter);

  qpb_spinor_ax(sum, constant_term, x);
  for(int sigma=0; sigma<number_of_terms; sigma++)
  {
    qpb_complex complex_numerator = {numerators[sigma], 0.};
    qpb_spinor_axpy(sum, complex_numerator, yMS[sigma], sum);
  }

  // // Testing results
  // print("\tTesting multi-sift solver results:\n");
  // qpb_double b_norm, yMS_norm = 0.;
  // qpb_spinor_xdotx(&b_norm, x);
  // // print("b_norm= %e\n ", b_norm);
  // for(int sigma=0; sigma<number_of_terms; sigma++)
  // {
  //   XdaggerX_plus_shift_op(y, yMS[sigma], shifts[sigma]);
  //   qpb_spinor_xmy(vector_diff, y, x);
  //   qpb_spinor_xdotx(&yMS_norm, vector_diff);
  //   print(" \t (||y[%d]-b||)/||b|| = %e\n", sigma, yMS_norm/b_norm);
  // }
  // print("\n");

  D_op(y, sum);

  qpb_mscongrad_finalize(number_of_terms);
}


void
qpb_overlap_kl(qpb_spinor_field y, qpb_spinor_field x, \
  enum qpb_kl_classes kl_class, int n, qpb_double epsilon, int max_iter)
{
  // Calculate the numerical terms of the partial fraction expression

  qpb_complex constant_term = {0., 0.};
  qpb_double *numerators, *shifts;
  numerators = qpb_alloc(sizeof(qpb_double)*n);
  shifts = qpb_alloc(sizeof(qpb_double)*n);

  constant_term.re = 1./(qpb_double) (2*n+1);
  for(int i=0; i<n; i++)
  {
    numerators[i] = 2*constant_term.re/pow(\
                              cos((M_PI/2.)*(2*i+1)/(qpb_double) (2*n+1)), 2);
    shifts[i] = pow(tan((M_PI/2.)*((2*i+1)/(qpb_double) (2*n+1))), 2);
  }

  // Modify the numerical terms using the scaling parameter
  constant_term.re *= 1/sqrt(scaling_factor);
  for(int i=0; i<n; i++)
  {
    numerators[i] *= sqrt(scaling_factor);
    shifts[i] *= scaling_factor;
  }

  /* Implementing (rho+overlap_mass/2)*x + (rho-overlap_mass/2)*g5(sign(X)) */
  qpb_double overlap_mass = ov_params.mass;
  qpb_double rho = ov_params.rho;

  qpb_spinor_field z = ov_temp_vecs[5];
  qpb_spinor_field w = ov_temp_vecs[6];

  qpb_gamma5_sign_function_partial_fractions(z, x, constant_term,\
                        numerators, shifts, n, epsilon, max_iter);

  qpb_complex a = {rho + 0.5*overlap_mass, 0.};
  qpb_complex b = {rho - 0.5*overlap_mass, 0.};

  qpb_spinor_axpby(y, a, x, b, z);

  free(numerators);
  free(shifts);

  return;
}


// ---------------------------- Single fraction ---------------------------- //

INLINE void
gamma5_single_fraction_numerator_op(qpb_spinor_field z, qpb_spinor_field x,\
                        unsigned long* single_fraction_coefficients, int n)
{
  qpb_spinor_field X_daggerX_to_power = ov_temp_vecs[7];
  qpb_spinor_field sum = ov_temp_vecs[8];

  // Initialize sum as 0
  qpb_spinor_field_set_zero(sum);
  // 0th-power term
  qpb_spinor_xeqy(X_daggerX_to_power, x);
  for(int i=0; i<n; i++)
  {
    qpb_complex complex_coefficient = {single_fraction_coefficients[i]\
                                              *pow(scaling_factor, n-i), 0.0};
    qpb_spinor_axpy(sum, complex_coefficient, X_daggerX_to_power, sum);
    XdaggerX_op(X_daggerX_to_power, X_daggerX_to_power);
  }
  // At the last term to the sum without any coefficient
  qpb_spinor_xpy(sum, X_daggerX_to_power, sum);

  // Calculate g5(sign(X))
  D_op(z, sum);

  return;
}


INLINE void
single_fraction_denominator_op(qpb_spinor_field sum, qpb_spinor_field x,\
                        unsigned long* single_fraction_coefficients, int n)
{
  qpb_spinor_field X_daggerX_to_power = ov_temp_vecs[9];
  qpb_complex reverse_complex_coefficient = {0.0, 0.0};

  // Initialize sum with the first term which is simply x
  // qpb_spinor_xeqy(sum, x);
  reverse_complex_coefficient.re = pow(scaling_factor, n + 1.0/2.0);
  qpb_spinor_ax(sum, reverse_complex_coefficient, x);
  // w will contain powers of x
  qpb_spinor_xeqy(X_daggerX_to_power, x);
  for(int i=0; i<n; i++)
  {
    XdaggerX_op(X_daggerX_to_power, X_daggerX_to_power);
    // Visit array elements in reverse order
    reverse_complex_coefficient.re = single_fraction_coefficients[n-i-1]\
                                      *pow(scaling_factor, (n-i-1) + 1.0/2.0);
    qpb_spinor_axpy(sum, reverse_complex_coefficient, X_daggerX_to_power, sum);
  }

  return;
}


// Function to calculate factorial
unsigned long factorial(int num)
{
  if (num == 0 || num == 1)
    return 1;
  else
    return num * factorial(num - 1);
}

// Function to calculate combinations (n choose k)
unsigned long calculateCombination(int n, int k)
{
    if (n < k)
      return 0;

    return factorial(n) / (factorial(k) * factorial(n - k));
}


int
qpb_congrad_single_fraction_denominator(qpb_spinor_field x, qpb_spinor_field b,\
                      unsigned long* single_fraction_coefficients, int n,\
                              qpb_double epsilon, int max_iter, int n_echo)
{
  qpb_spinor_field p = cg_temp_vecs[3];
  qpb_spinor_field r = cg_temp_vecs[4];
  qpb_spinor_field y = cg_temp_vecs[5];

  int n_reeval = 100;
  int iters = 0;

  qpb_double res_norm, b_norm;
  qpb_complex_double alpha = {1, 0}, omega = {1, 0};
  qpb_complex_double beta, gamma;
  qpb_spinor_xdotx(&b_norm, b);
  
  /* r = Denominator(x) */
  single_fraction_denominator_op(p, x, single_fraction_coefficients, n);
  /* r = b - Denominator(x) */
  qpb_spinor_xmy(r, b, p);
  qpb_spinor_xdotx(&gamma.re, r);
  gamma.im = 0;
  res_norm = gamma.re;
  qpb_spinor_xeqy(p, r);
  qpb_double t = qpb_stop_watch(0);
  for(iters=1; iters<max_iter; iters++)
  {
    if(res_norm / b_norm <= epsilon)
	    break;
    /* y = Denominator(p) */
    single_fraction_denominator_op(y, p, single_fraction_coefficients, n);

    /* omega = dot(p, Denominator(p)) */
    qpb_spinor_xdoty(&omega, p, y);

    /* alpha = dot(r, r)/omega */
    alpha = CDEV(gamma, omega);

    /* x <- x + alpha*p */
    qpb_spinor_axpy(x, alpha, p, x);
    if(iters % n_reeval == 0) 
    {
      A_op(y, x);
      qpb_spinor_xmy(r, b, y);
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
  single_fraction_denominator_op(y, x, single_fraction_coefficients, n);
  qpb_spinor_xmy(r, b, y);
  qpb_spinor_xdotx(&res_norm, r);

  if(iters==max_iter)
  {
    error(" !\n");
    error(" CG *did not* converge, after %d iterations\n", iters);
    error(" residual = %e, relative = %e, t = %g secs\n", res_norm,\
                                                        res_norm / b_norm, t);
    error(" !\n");
    return -1;
  }

  print(" \t After %d iters, CG converged, res = %e, relative = %e,\
                        t = %g secs\n", iters, res_norm, res_norm / b_norm, t);
  
  return iters;
}


// qpb_overlap_kl_single_fraction
// qpb_overlap_kl
void
qpb_overlap_kl_single_fraction(qpb_spinor_field y, qpb_spinor_field x,\
        enum qpb_kl_classes kl_class, int n, qpb_double epsilon, int max_iter)
{
  // Calculate the numerical coefficients for the single fraction expression
  unsigned long *single_fraction_coefficients;
  single_fraction_coefficients = qpb_alloc(sizeof(unsigned long)*n);

  for(int i=0; i<n; i++)
    single_fraction_coefficients[i] = calculateCombination(2*n+1, 2*i+1);

  // Calculate the denominator
	int n_echo = 100;
  qpb_spinor_field z = ov_temp_vecs[10];
  qpb_spinor_field_set_zero(z);

  qpb_congrad_single_fraction_denominator(z, x, single_fraction_coefficients,\
                                                 n, epsilon, max_iter, n_echo);

  // Calculate the numerator
  qpb_spinor_field w = ov_temp_vecs[11];
  gamma5_single_fraction_numerator_op(w, z, single_fraction_coefficients, n);

  /* Implementing (rho+overlap_mass/2)*x + (rho-overlap_mass/2)*g5(sign(X)) */
  qpb_double overlap_mass = ov_params.mass;
  qpb_double rho = ov_params.rho;

  qpb_complex a = {rho + 0.5*overlap_mass, 0.};
  qpb_complex b = {rho - 0.5*overlap_mass, 0.};

  qpb_spinor_axpby(y, a, x, b, w);

  free(single_fraction_coefficients);

  return;
}


// ---------------------------- Propagator ---------------------------- //

void
qpb_gamma5_overlap_kl(qpb_spinor_field y, qpb_spinor_field x, \
  enum qpb_kl_classes kl_class, int kl_iters, qpb_double epsilon, int max_iter)
{
  qpb_overlap_kl(y, x, kl_class, kl_iters, epsilon, max_iter);
  qpb_spinor_gamma5(y, y);

  return;
}


int
qpb_congrad_overlap_kl(qpb_spinor_field x, qpb_spinor_field b,\
  enum qpb_kl_classes kl_class, int kl_iters, qpb_double epsilon, int max_iter)
{
  qpb_spinor_field p = ov_temp_vecs[12];
  qpb_spinor_field r = ov_temp_vecs[13];
  qpb_spinor_field y = ov_temp_vecs[14];
  qpb_spinor_field w = ov_temp_vecs[15];
  qpb_spinor_field bprime = ov_temp_vecs[16];

  int n_reeval = 100;
  int n_echo = 100;
  int iters = 0;

  qpb_double res_norm, b_norm;
  qpb_complex_double alpha = {1, 0}, omega = {1, 0};
  qpb_complex_double beta, gamma;

  qpb_spinor_gamma5(w, b);
  qpb_gamma5_overlap_kl(bprime, w, kl_class, kl_iters, epsilon, max_iter);

  qpb_spinor_xdotx(&b_norm, bprime);

  /* r0 = bprime - A(x) */
  qpb_spinor_field_set_zero(p);
  // qpb_gamma5_overlap_kl(w, x, kl_class, kl_iters, epsilon, max_iter);
  // qpb_gamma5_overlap_kl(p, w, kl_class, kl_iters, epsilon, max_iter);
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
    qpb_gamma5_overlap_kl(w, p, kl_class, kl_iters, epsilon, max_iter);
    qpb_gamma5_overlap_kl(y, w, kl_class, kl_iters, epsilon, max_iter);

    /* omega = dot(p, A(p)) */
    qpb_spinor_xdoty(&omega, p, y);

    /* alpha = dot(r, r)/omega */
    alpha = CDEV(gamma, omega);

    /* x <- x + alpha*p */
    qpb_spinor_axpy(x, alpha, p, x);
    if(iters % n_reeval == 0) 
    {
      qpb_gamma5_overlap_kl(w, x, kl_class, kl_iters, epsilon, max_iter);
      qpb_gamma5_overlap_kl(y, w, kl_class, kl_iters, epsilon, max_iter);
      qpb_spinor_xmy(r, b, y);
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
  qpb_gamma5_overlap_kl(w, x, kl_class, kl_iters, epsilon, max_iter);
  qpb_gamma5_overlap_kl(y, w, kl_class, kl_iters, epsilon, max_iter);
  qpb_spinor_xmy(r, b, y);
  qpb_spinor_xdotx(&res_norm, r);

  if(iters==max_iter)
  {
    error(" !\n");
    error(" CG *did not* converge, after %d iterations\n", iters);
    error(" residual = %e, relative = %e, t = %g secs\n", res_norm,\
                                                        res_norm / b_norm, t);
    error(" !\n");
    return -1;
  }

  print(" \t After %d iters, CG converged, res = %e, relative = %e,\
                        t = %g secs\n", iters, res_norm, res_norm / b_norm, t);
  
  return iters;
}


INLINE void
propagator_single_fraction_denominator_op(qpb_spinor_field y, qpb_spinor_field x,\
                        unsigned long* single_fraction_coefficients, int n)
{
  qpb_spinor_field numerator = ov_temp_vecs[17];
  qpb_spinor_field denominator = ov_temp_vecs[18];

  qpb_double overlap_mass = ov_params.mass;
  qpb_double rho = ov_params.rho;

  qpb_complex a = {rho + 0.5*overlap_mass, 0.};
  qpb_complex b = {rho - 0.5*overlap_mass, 0.};

  gamma5_single_fraction_numerator_op(numerator, x, single_fraction_coefficients, n);
  single_fraction_denominator_op(denominator, x, single_fraction_coefficients, n);

  qpb_spinor_axpby(y, a, denominator, b, numerator);

  return;
}


int
qpb_congrad_propagator_single_fraction_denominator_op(qpb_spinor_field x, qpb_spinor_field b,\
                      unsigned long* single_fraction_coefficients, int n,\
                              qpb_double epsilon, int max_iter, int n_echo)
{  
  qpb_spinor_field p = cg_temp_vecs[6];
  qpb_spinor_field r = cg_temp_vecs[7];
  qpb_spinor_field y = cg_temp_vecs[8];
  qpb_spinor_field temp = cg_temp_vecs[9];

  int n_reeval = 100;
  int iters = 0;

  qpb_double res_norm, b_norm;
  qpb_complex_double alpha = {1, 0}, omega = {1, 0};
  qpb_complex_double beta, gamma;

  qpb_spinor_gamma5(b, b);
  propagator_single_fraction_denominator_op(y, b, single_fraction_coefficients, n);
  qpb_spinor_gamma5(b, y);
  qpb_spinor_xdotx(&b_norm, b);
  
  /* r = Denominator(x) */
  // propagator_single_fraction_denominator_op(p, x, single_fraction_coefficients, n);
  qpb_spinor_field_set_zero(p);

  /* r = b - Denominator(x) */
  qpb_spinor_xmy(r, b, p);
  qpb_spinor_xdotx(&gamma.re, r);
  gamma.im = 0;
  res_norm = gamma.re;
  qpb_spinor_xeqy(p, r);
  qpb_double t = qpb_stop_watch(0);
  for(iters=1; iters<max_iter; iters++)
  {
    if(res_norm / b_norm <= epsilon)
	    break;
    /* y = Denominator(p) */
    propagator_single_fraction_denominator_op(temp, p, single_fraction_coefficients, n);
    qpb_spinor_gamma5(y, temp);
    propagator_single_fraction_denominator_op(temp, y, single_fraction_coefficients, n);
    qpb_spinor_gamma5(y, temp);

    /* omega = dot(p, Denominator(p)) */
    qpb_spinor_xdoty(&omega, p, y);

    /* alpha = dot(r, r)/omega */
    alpha = CDEV(gamma, omega);

    /* x <- x + alpha*p */
    qpb_spinor_axpy(x, alpha, p, x);
    if(iters % n_reeval == 0) 
    {
      propagator_single_fraction_denominator_op(temp, x, single_fraction_coefficients, n);
      qpb_spinor_gamma5(y, temp);
      propagator_single_fraction_denominator_op(temp, y, single_fraction_coefficients, n);
      qpb_spinor_gamma5(y, temp);

      qpb_spinor_xmy(r, b, y);
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
  propagator_single_fraction_denominator_op(temp, x, single_fraction_coefficients, n);
  qpb_spinor_gamma5(y, temp);
  propagator_single_fraction_denominator_op(temp, y, single_fraction_coefficients, n);
  qpb_spinor_gamma5(y, temp);

  qpb_spinor_xmy(r, b, y);
  qpb_spinor_xdotx(&res_norm, r);

  if(iters==max_iter)
  {
    error(" !\n");
    error(" CG *did not* converge, after %d iterations\n", iters);
    error(" residual = %e, relative = %e, t = %g secs\n", res_norm,\
                                                        res_norm / b_norm, t);
    error(" !\n");
    return -1;
  }

  print(" \t After %d iters, CG converged, res = %e, relative = %e,\
                        t = %g secs\n", iters, res_norm, res_norm / b_norm, t);
  
  return iters;
}


int
qpb_inverse_overlap_kl_single_fraction(qpb_spinor_field sum,\
                qpb_spinor_field x, enum qpb_kl_classes kl_class, int n,\
                                            qpb_double epsilon, int max_iter)
{
  // Calculate the numerical coefficients for the single fraction expression
  unsigned long *single_fraction_coefficients;
  single_fraction_coefficients = qpb_alloc(sizeof(unsigned long)*n);

  for(int i=0; i<n; i++)
    single_fraction_coefficients[i] = calculateCombination(2*n+1, 2*i+1);

  // Calculate the denominator
	int n_echo = 100;
  qpb_spinor_field z = ov_temp_vecs[19];
  qpb_spinor_field_set_zero(z);

  int iters = qpb_congrad_propagator_single_fraction_denominator_op(z, x,\
                  single_fraction_coefficients, n, epsilon, max_iter, n_echo);


  // Calculate the numerator
  single_fraction_denominator_op(sum, z, single_fraction_coefficients, n);

  free(single_fraction_coefficients);

  return iters;
}
