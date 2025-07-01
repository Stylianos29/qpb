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


#define OVERLAP_NUMB_TEMP_VECS 11
#define CG_NUMB_TEMP_VECS 13

static qpb_spinor_field ov_temp_vecs[OVERLAP_NUMB_TEMP_VECS];
static qpb_spinor_field cg_temp_vecs[CG_NUMB_TEMP_VECS];

static qpb_overlap_params ov_params;

static int KL_diagonal_order;
static qpb_double inner_solver_precision;
static int maximum_inner_solver_iterations;
static qpb_double scaling_factor;

unsigned long *sum_coefficients_sfrac;


// HELPER FUNCTIONS

unsigned long factorial(int num)
{
  /* Function to calculate factorial */
  if (num == 0 || num == 1)
    return 1;
  else
    return num * factorial(num - 1);
}


unsigned long calculate_combination(int n, int k)
{
  /* Function to calculate combinations (n choose k) */
  if (n < k)
    return 0;

  return factorial(n) / (factorial(k) * factorial(n - k));
}


void
qpb_overlap_kl_sfrac_init(void * gauge, qpb_clover_term clover, \
          enum qpb_kl_classes kl_class, int kl_iters, qpb_double rho, \
          qpb_double c_sw, qpb_double mass, qpb_double mu, \
          qpb_double inner_epsilon, int max_inner_iter)
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

    KL_diagonal_order = kl_iters;
    inner_solver_precision = inner_epsilon;
    maximum_inner_solver_iterations = max_inner_iter;
    scaling_factor = mu;

    // Calculate the numerical coefficients for the numerator of the single
    // fraction expression. The denominator coefficients mirror the numerator's
    // but are arranged in reverse order.
    sum_coefficients_sfrac = qpb_alloc(sizeof(unsigned long)\
                                                        *KL_diagonal_order+1);
    
    for(int i=0; i<KL_diagonal_order+1; i++)
    {
      sum_coefficients_sfrac[i] = calculate_combination(2*KL_diagonal_order+1, \
                                                                        2*i+1);
      // print("sum_coefficients_sfrac[%d] = %d\n", i, sum_coefficients_sfrac[i]);
                                                                    
    }
  }
	
  return;
}


void
qpb_overlap_kl_sfrac_finalize()
{
  qpb_comm_halo_spinor_field_finalize();
  for(int i=0; i<OVERLAP_NUMB_TEMP_VECS; i++)
    qpb_spinor_field_finalize(ov_temp_vecs[i]);
  
  for(int i=0; i<CG_NUMB_TEMP_VECS; i++)
    qpb_spinor_field_finalize(cg_temp_vecs[i]);

  if(which_dslash_op == QPB_DSLASH_STANDARD)
    qpb_gauge_field_finalize(*(qpb_gauge_field *)ov_params.gauge_ptr);
  
  ov_params.initialized = 0;
  
  free(sum_coefficients_sfrac);
  
  return;
}


INLINE void
D_op(qpb_spinor_field y, qpb_spinor_field x)
{
  /* Implements y = (D - rho*I)(x) */

  void *dslash_args[4];

  dslash_args[0] = ov_params.gauge_ptr;
  dslash_args[1] = &ov_params.m_bare;
  dslash_args[2] = &ov_params.clover;
  dslash_args[3] = &ov_params.c_sw;
  
  ov_params.dslash_op(y, x, dslash_args);
  
  return;
}


INLINE void
X2_op(qpb_spinor_field y, qpb_spinor_field x)
{
  /* Implements y = X^2(x) = γ5(D - rho*I)(γ5(D - rho*I))(x) */

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


INLINE void
pnn_X2_op(qpb_spinor_field y, qpb_spinor_field x)
{
  /* Implements: y = p_{nn}(X^2(x)), with:
  - X(x) = γ5(D(x) - ρ*x) = (D - rho*I)(x), 
  - p_{nn}(X^2(x)) = sum_{k=0}^n C(2*n+1, 2*k+1) (X^2)^k(x) */

  qpb_spinor_field X2_to_power = ov_temp_vecs[1];
  qpb_spinor_field sum = ov_temp_vecs[2];

  // Initialize sum to zero
  qpb_spinor_field_set_zero(sum);
  // 0th-power term
  qpb_spinor_xeqy(X2_to_power, x);
  for(int i=0; i<KL_diagonal_order; i++)
  {
    qpb_complex complex_coefficient = {sum_coefficients_sfrac[i]\
                              *pow(scaling_factor, KL_diagonal_order-i), 0.0};
    qpb_spinor_axpy(sum, complex_coefficient, X2_to_power, sum);
    X2_op(X2_to_power, X2_to_power);
  }
  // Add the last term to the sum without any coefficient
  qpb_spinor_xpy(y, X2_to_power, sum);

  return;
}


INLINE void
qnn_X2_op(qpb_spinor_field y, qpb_spinor_field x)
{
  /* Implements: y = q_{nn}(X^2(x)), with:
  - X(x) = γ5(D(x) - ρ*x) = (D - rho*I)(x), 
  - q_{nn}(X^2(x)) = sum_{k=0}^n C(2*n+1, 2*k) (X^2)^k(x) */

  qpb_spinor_field X2_to_power = ov_temp_vecs[3];
  qpb_spinor_field sum = ov_temp_vecs[4];
  qpb_complex complex_coefficient = {0.0, 0.0};

  // Initialize sum with the first term which is simply x
  complex_coefficient.re = 1.0;
  if (scaling_factor != 1.0)
    complex_coefficient.re *= pow(scaling_factor, KL_diagonal_order + 0.5);

  // print("complex_coefficient.re = %g\n", complex_coefficient.re);
  qpb_spinor_ax(sum, complex_coefficient, x);
  qpb_spinor_xeqy(X2_to_power, x);
  for(int i=0; i<KL_diagonal_order; i++)
  {
    X2_op(X2_to_power, X2_to_power);
    // Visit array elements in reverse order
    complex_coefficient.re = sum_coefficients_sfrac[KL_diagonal_order-i-1];
    if (scaling_factor != 1.0)
      complex_coefficient.re *= pow(scaling_factor, (KL_diagonal_order-i-1) + 0.5);
    qpb_spinor_axpy(sum, complex_coefficient, X2_to_power, sum);
  }
  qpb_spinor_xeqy(y, sum);

  return;
}


int
qpb_congrad_qnn_X2(qpb_spinor_field x, qpb_spinor_field b)
{
  qpb_spinor_field p = cg_temp_vecs[0];
  qpb_spinor_field r = cg_temp_vecs[1];
  qpb_spinor_field y = cg_temp_vecs[2];

  int n_reeval = 100;
  int n_echo = 100;
  int iters = 0;

  qpb_double res_norm, b_norm;
  qpb_complex_double alpha = {1, 0}, omega = {1, 0};
  qpb_complex_double beta, gamma;

  qpb_spinor_xdotx(&b_norm, b);

  qpb_spinor_field_set_zero(x);

  /* r0 = b - A(x) */
  // qnn_X2_op(p, x);
  // qpb_spinor_xmy(r, b, p);

  /* Or r0 = bprime for short since x0 = 0 */
  qpb_spinor_xeqy(r, b);

  qpb_spinor_xdotx(&gamma.re, r);
  gamma.im = 0;
  res_norm = gamma.re;
  /* p = r0 */
  qpb_spinor_xeqy(p, r);

  qpb_double t = qpb_stop_watch(0);
  for(iters=1; iters<maximum_inner_solver_iterations; iters++)
  {
    if(res_norm / b_norm <= inner_solver_precision)
    {
      // print("Inner CG stopped at relative residual: %e\n", res_norm / b_norm);
      break;
    }

    /* y = A(p) */
    qnn_X2_op(y, p);

    /* omega = dot(p, A(p)) */
    qpb_spinor_xdoty(&omega, p, y);

    /* alpha = dot(r, r)/omega */
    alpha = CDEV(gamma, omega);

    /* x <- x + alpha*p */
    qpb_spinor_axpy(x, alpha, p, x);

    if(iters % n_reeval == 0) 
    {
      qnn_X2_op(y, x);
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
  qnn_X2_op(y, x);
  qpb_spinor_xmy(r, b, y);
  qpb_spinor_xdotx(&res_norm, r);

  if(iters==maximum_inner_solver_iterations)
  {
    error(" !\n");
    error(" q_nn CG *did not* converge, after %d iterations\n", iters);
    error(" residual = %e, relative = %e, t = %g sec\n", res_norm,\
                                                        res_norm / b_norm, t);
    error(" !\n");
    return -1;
  }

  // print(" \tAfter %d iters, q_nn CG converged, res = %e, "
  //         "relative = %e, t = %g sec\n",
  //          iters, res_norm, res_norm / b_norm, t);
  
  return iters;
}


void
qpb_gamma5_sign_function_of_X_sfrac(qpb_spinor_field y, qpb_spinor_field x)
{
  /* Implements: y = γ5(sign(X(x))) = γ5(X(p_{nn}(X^2)(q_{nn}(X^2)^-1)))(x),
  with:
  - X(x) = γ5(D(x) - ρ*x) = γ5(D - rho*I)(x), 
  - p_{nn}(X^2(x)) = sum_{k=0}^n C(2*n+1, 2*k+1) (X^2)^k(x)
  - q_{nn}(X^2(x)) = sum_{k=0}^n C(2*n+1, 2*k) (X^2)^k(x)
   */

  qpb_spinor_field z = ov_temp_vecs[5];
  qpb_spinor_field w = ov_temp_vecs[6];

  qpb_congrad_qnn_X2(w, x);
  pnn_X2_op(z, w);
  // y = γ5(X(z)) = (D - rho*I)(z)
  D_op(y, z);

  return;
}


void
qpb_overlap_kl_sfrac(qpb_spinor_field y, qpb_spinor_field x)
{
  /* Implements:
        Dov,m(x) = (rho+overlap_mass/2)*x+ (rho-overlap_mass/2)*γ5(sign(X))
  */
  
  qpb_spinor_field z = ov_temp_vecs[7];

  qpb_double overlap_mass = ov_params.mass; // Overlap operator Dov,m mass
  qpb_double rho = ov_params.rho;

  qpb_complex a = {rho + 0.5*overlap_mass, 0.};
  qpb_complex b = {rho - 0.5*overlap_mass, 0.};

  qpb_gamma5_sign_function_of_X_sfrac(z, x);

  qpb_spinor_axpby(y, a, x, b, z);

  return;
}


void
qpb_gamma5_overlap_kl_sfrac(qpb_spinor_field y, qpb_spinor_field x)
{
  qpb_overlap_kl_sfrac(y, x);
  qpb_spinor_gamma5(y, y);

  return;
}


int
qpb_congrad_overlap_kl_sfrac(qpb_spinor_field x, qpb_spinor_field b, \
                                        qpb_double CG_epsilon, int CG_max_iter)
{
  qpb_spinor_field p = cg_temp_vecs[3];
  qpb_spinor_field r = cg_temp_vecs[4];
  qpb_spinor_field y = cg_temp_vecs[5];
  qpb_spinor_field w = cg_temp_vecs[6];
  qpb_spinor_field bprime = cg_temp_vecs[7];

  int n_reeval = 100;
  int n_echo = 10;
  int iters = 0;

  qpb_double res_norm, b_norm;
  qpb_complex_double alpha = {1, 0}, omega = {1, 0};
  qpb_complex_double beta, gamma;

  // bprime = γ5(D_ov(γ5(b)))
  qpb_spinor_gamma5(w, b);
  qpb_gamma5_overlap_kl_sfrac(bprime, w);

  qpb_spinor_xdotx(&b_norm, bprime);

  qpb_spinor_field_set_zero(x);

  /* r0 = bprime - A(x) */
  // qpb_gamma5_overlap_kl_sfrac(w, x);
  // qpb_gamma5_overlap_kl_sfrac(p, w);
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
    if(res_norm / b_norm <= CG_epsilon)
      {
        // print("CG stopped at relative residual: %e\n", res_norm / b_norm);
        break;
      }

    /* y = A(p) */
    qpb_gamma5_overlap_kl_sfrac(w, p);
    qpb_gamma5_overlap_kl_sfrac(y, w);

    /* omega = dot(p, A(p)) */
    qpb_spinor_xdoty(&omega, p, y);

    /* alpha = dot(r, r)/omega */
    alpha = CDEV(gamma, omega);

    /* x <- x + alpha*p */
    qpb_spinor_axpy(x, alpha, p, x);

    if(iters % n_reeval == 0) 
    {
      qpb_gamma5_overlap_kl_sfrac(w, x);
      qpb_gamma5_overlap_kl_sfrac(y, w);
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

  qpb_gamma5_overlap_kl_sfrac(w, x);
  qpb_gamma5_overlap_kl_sfrac(y, w);
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

  // print(" \tAfter %d iters, CG converged, res = %e, "
  //           "relative = %e, t = %g sec\n",
  //            iters, res_norm, res_norm / b_norm, t);

  return iters;
}

/* --------------------- IMPLEMENTING MULTIPLY-UP TRICK --------------------- */

void
aQnn_gamma5_plus_bX_Pnn_op(qpb_spinor_field y, qpb_spinor_field x)
{

  qpb_spinor_field z = ov_temp_vecs[8];
  qpb_spinor_field w = ov_temp_vecs[9];

  qpb_double overlap_mass = ov_params.mass; // Overlap operator Dov,m mass
  qpb_double rho = ov_params.rho;

  qpb_complex a = {rho + 0.5*overlap_mass, 0.};
  qpb_complex b = {rho - 0.5*overlap_mass, 0.};

  pnn_X2_op(y, x);
  D_op(z, y);
  qpb_spinor_gamma5(z, z);

  qpb_spinor_gamma5(y, x);
  qnn_X2_op(w, y);

  qpb_spinor_axpby(y, a, w, b, z);

  return;
}


void
conjugate_aQnn_gamma5_plus_bX_Pnn_op(qpb_spinor_field y, qpb_spinor_field x)
{

  qpb_spinor_field z = ov_temp_vecs[8];
  qpb_spinor_field w = ov_temp_vecs[9];

  qpb_double overlap_mass = ov_params.mass; // Overlap operator Dov,m mass
  qpb_double rho = ov_params.rho;

  qpb_complex a = {rho + 0.5*overlap_mass, 0.};
  qpb_complex b = {rho - 0.5*overlap_mass, 0.};

  pnn_X2_op(y, x);
  D_op(z, y);
  qpb_spinor_gamma5(z, z);

  qnn_X2_op(y, x);
  qpb_spinor_gamma5(w, y);

  qpb_spinor_axpby(y, a, w, b, z);

  return;
}


int
qpb_congrad_aqnn_plus_bDpnn(qpb_spinor_field x, \
                    qpb_spinor_field b, qpb_double CG_epsilon, int CG_max_iter)
{
  qpb_spinor_field p = cg_temp_vecs[8];
  qpb_spinor_field r = cg_temp_vecs[9];
  qpb_spinor_field y = cg_temp_vecs[10];
  qpb_spinor_field w = cg_temp_vecs[11];
  qpb_spinor_field bprime = cg_temp_vecs[12];

  int n_reeval = 100;
  int n_echo = 100;
  int iters = 0;

  qpb_double res_norm, b_norm;
  qpb_complex_double alpha = {1, 0}, omega = {1, 0};
  qpb_complex_double beta, gamma;

  // bprime = Dov^+ b
  // qpb_spinor_gamma5(w, b);
  conjugate_aQnn_gamma5_plus_bX_Pnn_op(bprime, b);

  qpb_spinor_xdotx(&b_norm, bprime);
  
  qpb_spinor_field_set_zero(x);

  /* r0 = bprime - A(x) */
  // aQnn_gamma5_plus_bX_Pnn_op(w, x);
  // conjugate_aQnn_gamma5_plus_bX_Pnn_op(p, w);
  // qpb_spinor_xmy(r, bprime, p);

  /* Or r0 = bprime for short since x0=0 */
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
      // print("CG stopped at relative residual: %e\n", res_norm / b_norm);
      break;
    }

    /* y = A(p) */
    aQnn_gamma5_plus_bX_Pnn_op(w, p);
    conjugate_aQnn_gamma5_plus_bX_Pnn_op(y, w);

    /* omega = dot(p, A(p)) */
    qpb_spinor_xdoty(&omega, p, y);

    /* alpha = dot(r, r)/omega */
    alpha = CDEV(gamma, omega);

    /* x <- x + alpha*p */
    qpb_spinor_axpy(x, alpha, p, x);

    if(iters % n_reeval == 0) 
    {
      aQnn_gamma5_plus_bX_Pnn_op(w, x);
      conjugate_aQnn_gamma5_plus_bX_Pnn_op(y, w);
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

  aQnn_gamma5_plus_bX_Pnn_op(w, x);
  conjugate_aQnn_gamma5_plus_bX_Pnn_op(y, w);
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

  printf(" \tAfter %d iters, CG converged, res = %e, "
          "relative = %e, t = %g sec\n",
           iters, res_norm, res_norm / b_norm, t);

  return iters;
}


int
qpb_congrad_overlap_kl_sfrac_multiply_up(qpb_spinor_field x, \
                    qpb_spinor_field b, qpb_double CG_epsilon, int CG_max_iter)
{
  qpb_spinor_field b_prime = ov_temp_vecs[10];

  qpb_spinor_gamma5(b, b);
  qnn_X2_op(b_prime, b);

  qpb_congrad_aqnn_plus_bDpnn(x, b_prime, CG_epsilon, CG_max_iter);
}
