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

    /* Calculate the numerical terms of the partial fraction expansion */
    shifts = qpb_alloc(sizeof(qpb_double)*KL_diagonal_order);
    numerators = qpb_alloc(sizeof(qpb_double)*KL_diagonal_order);

    constant_term = 1.0/((qpb_double) (2*KL_diagonal_order+1));

    for(int i=0; i<KL_diagonal_order; i++)
    {
      qpb_double trig_arg = M_PI*(i+0.5)*constant_term;
      shifts[i] = pow(tan(trig_arg), 2);
      numerators[i] = 2*constant_term/powl(cos(trig_arg), 2);
      // print("numerator[%d] = %.25f, shift[%d] = %.25f\n", i, numerators[i], \
                                                              i, shifts[i]);
    }

    // Modify the numerical constants of the partial fraction expansions using
    // the scaling parameter
    if (scaling_factor != 1.0)
    {
      constant_term *= 1/sqrt(scaling_factor);
      for(int i=0; i<KL_diagonal_order; i++)
      {
        numerators[i] *= sqrt(scaling_factor);
        shifts[i] *= scaling_factor;
      }
    }

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
qpb_gamma5_sign_function_of_X_pfrac(qpb_spinor_field y, qpb_spinor_field x)
{
  /* Implements: γ5(sign(X(x))) = γ5(X(c_0 + Sum_{i=1}^{n} c_i/(X^2+σ_i) )),
      with X(x) = γ5(D(x) - ρ*x) . */

  qpb_spinor_field sum = ov_temp_vecs[0];

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

  // Initialize sum with the constant term
  qpb_spinor_ax(sum, (qpb_complex) {constant_term, 0.}, x);
  // And then add the rest of the partial fraction terms
  for(int sigma=0; sigma<KL_diagonal_order; sigma++)
    qpb_spinor_axpy(sum, (qpb_complex) {numerators[sigma], 0.}, yMS[sigma], sum);

  D_op(y, sum);

  return;
}


void
qpb_overlap_kl_pfrac(qpb_spinor_field y, qpb_spinor_field x)
{
  /* Implements:
        Dov,m(x) = (rho+overlap_mass/2)*x+ (rho-overlap_mass/2)*g5(sign(X))
  */
  
  qpb_spinor_field z = ov_temp_vecs[1];

  qpb_double overlap_mass = ov_params.mass; // Overlap operator Dov,m mass
  qpb_double rho = ov_params.rho;

  qpb_complex a = {rho + 0.5*overlap_mass, 0.};
  qpb_complex b = {rho - 0.5*overlap_mass, 0.};

  qpb_gamma5_sign_function_of_X_pfrac(z, x);

  qpb_spinor_axpby(y, a, x, b, z);

  return;
}


void
qpb_gamma5_overlap_kl_pfrac(qpb_spinor_field y, qpb_spinor_field x)
{
  qpb_overlap_kl_pfrac(y, x);
  qpb_spinor_gamma5(y, y);

  return;
}


int
qpb_congrad_overlap_kl_pfrac(qpb_spinor_field x, qpb_spinor_field b, \
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
  qpb_gamma5_overlap_kl_pfrac(bprime, w);

  qpb_spinor_xdotx(&b_norm, bprime);

  qpb_spinor_field_set_zero(x);

  /* r0 = bprime - A(x) */
  // qpb_gamma5_overlap_kl_pfrac(w, x);
  // qpb_gamma5_overlap_kl_pfrac(p, w);
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
    qpb_gamma5_overlap_kl_pfrac(w, p);
    qpb_gamma5_overlap_kl_pfrac(y, w);

    /* omega = dot(p, A(p)) */
    qpb_spinor_xdoty(&omega, p, y);

    /* alpha = dot(r, r)/omega */
    alpha = CDEV(gamma, omega);

    /* x <- x + alpha*p */
    qpb_spinor_axpy(x, alpha, p, x);

    if(iters % n_reeval == 0) 
    {
      qpb_gamma5_overlap_kl_pfrac(w, x);
      qpb_gamma5_overlap_kl_pfrac(y, w);
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

  qpb_gamma5_overlap_kl_pfrac(w, x);
  qpb_gamma5_overlap_kl_pfrac(y, w);
  qpb_spinor_xmy(r, bprime, y);
  qpb_spinor_xdotx(&res_norm, r);

  if(iters==CG_max_iter)
  {
    error(" !\n");
    error(" CG *did not* converge, after %d iterations\n", iters);
    error(" residual = %e, relative = %e, t = %g secs\n", res_norm, \
                                                      res_norm / b_norm, t);
    error(" !\n");
    return -1;
  }

  print(" \tAfter %d iters, CG converged, res = %e, relative = %e, "
        "t = %g secs\n",
         iters, res_norm, res_norm / b_norm, t);
  
  return iters;
}
