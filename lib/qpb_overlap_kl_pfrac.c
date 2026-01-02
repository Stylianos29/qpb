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


#define OVERLAP_NUMB_TEMP_VECS 9
#define MSCG_NUMB_TEMP_VECS 20


static qpb_spinor_field ov_temp_vecs[OVERLAP_NUMB_TEMP_VECS];
static qpb_spinor_field mscg_temp_vecs[MSCG_NUMB_TEMP_VECS];

static qpb_overlap_params ov_params;

static int KL_diagonal_order;
static qpb_double MS_solver_precision;
static int MS_maximum_solver_iterations;

static qpb_double rho_plus;
static qpb_double rho_minus;

static qpb_double constant_term;
static qpb_double *c;

static qpb_double left_numerator;
static qpb_double left_denominator;
static qpb_double right_numerator;
static qpb_double right_denominator;


void
qpb_overlap_kl_pfrac_init(void * gauge, qpb_clover_term clover, \
          enum qpb_kl_classes kl_class, int kl_iters, qpb_double rho, \
          qpb_double c_sw, qpb_double mass, qpb_double scaling_factor, \
          qpb_double ms_epsilon, int ms_max_iter)
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

    rho_plus = rho + 0.5*ov_params.mass;
    rho_minus = rho - 0.5*ov_params.mass;
    
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

    c = qpb_alloc(sizeof(qpb_double)*2*KL_diagonal_order);

    /* Calculate the numerical constants of the product form expansion
    of the sign function */
    constant_term = 1.0/((qpb_double) (2*KL_diagonal_order+1));
    for(int m=0; m<2*KL_diagonal_order; m++)
    {
      qpb_double trig_arg = (m+1)*constant_term*0.5*M_PI;
      c[m] = pow(tan(trig_arg), 2);
      // print("c[%d] = %.25f\n", c[m], m);
    }

    left_numerator = c[2];
    left_denominator = c[3];
    right_numerator = c[1];
    right_denominator = c[0];

    // TEMPORARY: Only initialize half the vectors for MSCG
    qpb_mscongrad_init(KL_diagonal_order);

  }
	
  return;
}


void
qpb_overlap_kl_pfrac_finalize()
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
  
  return;
}


INLINE void
X_op(qpb_spinor_field y, qpb_spinor_field x)
{
  /* X ≡ γ5 (a*D - ρ) */

  void *dslash_args[4];

  dslash_args[0] = ov_params.gauge_ptr;
  dslash_args[1] = &ov_params.m_bare;
  dslash_args[2] = &ov_params.clover;
  dslash_args[3] = &ov_params.c_sw;

  ov_params.g5_dslash_op(y, x, dslash_args);

  return;
}


void
qpb_first_degree_rational(qpb_spinor_field y, qpb_spinor_field x, qpb_double a, qpb_double b)
{
  /* Implements the partial fraction decomposition of a rational:
    (X^2 +a)/(X^2 + b) = 1 + (a-b)/(X^2 + b) .
  */

  qpb_double *shift = &b;

  qpb_spinor_field yMS[KL_diagonal_order/2];
  for(int sigma=0; sigma<KL_diagonal_order/2; sigma++)
  {
    yMS[sigma] = mscg_temp_vecs[sigma];
    // It needs to re-initialized to 0 with every call of the function
    qpb_spinor_field_set_zero(yMS[sigma]);
  }

  qpb_double kernel_mass = ov_params.m_bare; // Kernel operator bare mass
  qpb_double kernel_kappa = 1./(2*kernel_mass+8.);

  qpb_mscongrad(yMS, x, ov_params.gauge_ptr, ov_params.clover, kernel_kappa, \
    KL_diagonal_order/2, shift, ov_params.c_sw, MS_solver_precision, \
    MS_maximum_solver_iterations);

  qpb_spinor_axpy(y, (qpb_complex) {a-b, 0.}, yMS[0], x);

  return;
}


void
qpb_overlap_kl_pfrac_multiply_down(qpb_spinor_field y, qpb_spinor_field x)
{
  /* Implements: 
      
    (X^2 + c3)/(X^2 + c4) γ5 Dov,m(x) 
        = ρ+ (X^2 + c3)/(X^2 + c4) γ5 x + ρ- Χ (X^2 + c2)/(X^2 + c1)

    with ρ+ = ρ + overlap_mass/2 and ρ- = ρ - overlap_mass/2.  
  */

  qpb_spinor_field z = ov_temp_vecs[0];
  qpb_spinor_field w = ov_temp_vecs[1];

  qpb_spinor_gamma5(y, x);
  qpb_first_degree_rational(z, y, left_numerator, left_denominator);

  qpb_first_degree_rational(y, x, right_numerator, right_denominator);
  X_op(w, y);

  qpb_spinor_axpby(y, (qpb_complex) {rho_plus, 0.}, z, (qpb_complex) {rho_minus*constant_term, 0.}, w);

  return;
}


void
qpb_conjugate_overlap_kl_pfrac_multiply_down(qpb_spinor_field y, qpb_spinor_field x)
{
  /* Implements: 
      
      ρ+ γ5 (X^2 + c3)/(X^2 + c4) x + ρ- Χ (X^2 + c2)/(X^2 + c1)

    with ρ+ = ρ + overlap_mass/2 and ρ- = ρ - overlap_mass/2.  
  */

  qpb_spinor_field z = ov_temp_vecs[2];
  qpb_spinor_field w = ov_temp_vecs[3];

  qpb_double shifts_array[2] = {c[0], c[3]};
  qpb_double *shifts = shifts_array;

  // Initialize all MSCG temp vectors
  qpb_spinor_field yMS[KL_diagonal_order];
  for(int sigma=0; sigma<KL_diagonal_order; sigma++)
  {
    yMS[sigma] = mscg_temp_vecs[sigma];
    // It needs to re-initialized to 0 with every call of the function
    qpb_spinor_field_set_zero(yMS[sigma]);
  }

  qpb_double kernel_mass = ov_params.m_bare; // Kernel operator bare mass
  qpb_double kernel_kappa = 1./(2*kernel_mass+8.);

  qpb_mscongrad(yMS, x, ov_params.gauge_ptr, ov_params.clover, kernel_kappa, \
    KL_diagonal_order, shifts, ov_params.c_sw, MS_solver_precision, \
    MS_maximum_solver_iterations);

  qpb_spinor_axpy(y, (qpb_complex) {left_numerator - left_denominator, 0.}, yMS[1], x);
  qpb_spinor_gamma5(z, y);

  qpb_spinor_axpy(y, (qpb_complex) {right_numerator - right_denominator, 0.}, yMS[0], x);
  X_op(w, y);

  qpb_spinor_axpby(y, (qpb_complex) {rho_plus, 0.}, z, (qpb_complex) {rho_minus*constant_term, 0.}, w);

  return;
}


int
qpb_congrad_overlap_kl_pfrac(qpb_spinor_field x, qpb_spinor_field b, \
                                        qpb_double CG_epsilon, int CG_max_iter)
{
  qpb_spinor_field p = ov_temp_vecs[4];
  qpb_spinor_field r = ov_temp_vecs[5];
  qpb_spinor_field y = ov_temp_vecs[6];
  qpb_spinor_field w = ov_temp_vecs[7];
  qpb_spinor_field bprime = ov_temp_vecs[8];

  int n_reeval = 100;
  int n_echo = 100;
  int iters = 0;

  
  qpb_double res_norm, b_norm;
  qpb_complex_double alpha = {1, 0}, omega = {1, 0};
  qpb_complex_double beta, gamma;
  
  qpb_spinor_gamma5(b, b);
  qpb_first_degree_rational(w, b, left_numerator, left_denominator);

  qpb_conjugate_overlap_kl_pfrac_multiply_down(bprime, w);

  qpb_spinor_xdotx(&b_norm, bprime);

  qpb_spinor_field_set_zero(x);

  /* r0 = bprime - A(x) */
  // qpb_overlap_kl_pfrac_multiply_down(w, x);
  // qpb_conjugate_overlap_kl_pfrac_multiply_down(p, w);
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
    qpb_overlap_kl_pfrac_multiply_down(w, p);
    qpb_conjugate_overlap_kl_pfrac_multiply_down(y, w);

    /* omega = dot(p, A(p)) */
    qpb_spinor_xdoty(&omega, p, y);

    /* alpha = dot(r, r)/omega */
    alpha = CDEV(gamma, omega);

    /* x <- x + alpha*p */
    qpb_spinor_axpy(x, alpha, p, x);

    if(iters % n_reeval == 0) 
    {
      qpb_overlap_kl_pfrac_multiply_down(w, x);
      qpb_conjugate_overlap_kl_pfrac_multiply_down(y, w);
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

  qpb_overlap_kl_pfrac_multiply_down(w, x);
  qpb_conjugate_overlap_kl_pfrac_multiply_down(y, w);
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
