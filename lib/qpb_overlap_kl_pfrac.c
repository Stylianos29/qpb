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


#define OVERLAP_NUMB_TEMP_VECS 16
#define MSCG_NUMB_TEMP_VECS 20


static qpb_spinor_field ov_temp_vecs[OVERLAP_NUMB_TEMP_VECS];
static qpb_spinor_field mscg_temp_vecs[MSCG_NUMB_TEMP_VECS];

static qpb_overlap_params ov_params;

static int KL_diagonal_order;
static qpb_double MS_solver_precision;
static int MS_maximum_solver_iterations;

static qpb_double prec_CG_epsilon;
static int prec_CG_max_iter;

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

    prec_CG_epsilon = 1e-4;
    prec_CG_max_iter = 10000;

    print(" Preconditioner solver epsilon = %e\n", prec_CG_epsilon);
    print(" Preconditioner solver max iterations = %d\n", prec_CG_max_iter);

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

    /* Apply scaling parameter to the partial fraction coefficients */
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


INLINE void
M_op(qpb_spinor_field y, qpb_spinor_field x)
{

  // /* Implement M =   */ 

  qpb_spinor_field w = ov_temp_vecs[0];

  qpb_double overlap_mass = ov_params.mass;
  qpb_double rho = ov_params.rho;

  qpb_complex preconditioner_mass = {rho + overlap_mass, 0.};
  
  D_op(w, x);

  qpb_spinor_axpy(y, preconditioner_mass, x, w);
  
  return;
}


INLINE void
M_conj_op(qpb_spinor_field y, qpb_spinor_field x)
{

  qpb_spinor_field z = ov_temp_vecs[1];

  qpb_spinor_gamma5(y, x);
  M_op(z, y);
  qpb_spinor_gamma5(y, z);

  return;
}


void
qpb_gamma5_sign_function_of_X_pfrac(qpb_spinor_field y, qpb_spinor_field x)
{
  /* Implements: γ5(sign(X(x))) = γ5(X(c_0 + Sum_{i=1}^{n} c_i/(X^2+σ_i) )),
      with X(x) = γ5(D(x) - ρ*x) . */

  qpb_spinor_field sum = ov_temp_vecs[2];

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
        Dov,m(x) = (rho+overlap_mass/2)*x + (rho-overlap_mass/2)*g5(sign(X))
  */
  
  qpb_spinor_field z = ov_temp_vecs[3];

  qpb_double overlap_mass = ov_params.mass;
  qpb_double rho = ov_params.rho;

  qpb_complex rho_plus = {rho + 0.5*overlap_mass, 0.};
  qpb_complex rho_minus = {rho - 0.5*overlap_mass, 0.};

  qpb_gamma5_sign_function_of_X_pfrac(z, x);

  qpb_spinor_axpby(y, rho_plus, x, rho_minus, z);

  return;
}


void
qpb_conjugate_overlap_kl_pfrac(qpb_spinor_field y, qpb_spinor_field x)
{
  qpb_spinor_field z = ov_temp_vecs[4];

  qpb_spinor_gamma5(y, x);
  qpb_overlap_kl_pfrac(z, y);
  qpb_spinor_gamma5(y, z);

  return;
}


int
qpb_preconditioner_CG(qpb_spinor_field s, qpb_spinor_field b)
{
  qpb_spinor_field r = ov_temp_vecs[5];
  qpb_spinor_field p = ov_temp_vecs[6];
  qpb_spinor_field w = ov_temp_vecs[7];
  qpb_spinor_field y = ov_temp_vecs[8];

  int iters = 0;

  /* All scalars are real -- D†D is HPD */
  qpb_double b_norm, res_norm, new_res_norm, omega, alpha, beta;

  /* ||b||^2 */
  qpb_spinor_xdotx(&b_norm, b);

  /* s0 = 0,  r0 = b */
  qpb_spinor_field_set_zero(s);
  qpb_spinor_xeqy(r, b);

  /* gamma_0 = ||r0||^2 */
  qpb_spinor_xdotx(&res_norm, r);

  /* p0 = r0 */
  qpb_spinor_xeqy(p, r);

  for(iters = 1; iters < prec_CG_max_iter; iters++)
  {
    /* Stopping on relative residual of the HPD system:
       ||r||^2 / ||b||^2 <= prec_CG_epsilon^2   (squaring avoids a sqrt) */
    if(res_norm / b_norm <= prec_CG_epsilon * prec_CG_epsilon)
      break;

    /* Apply D†D to p: w = D·p,  y = D†·w = D†D·p */
    M_op(w, p);
    M_conj_op(y, w);

    /* omega = p†·D†D·p = ||D·p||^2 = ||w||^2  (real, positive) */
    qpb_spinor_xdotx(&omega, w);

    /* alpha = ||r||^2 / (p†·D†D·p)  (real, positive) */
    alpha = res_norm / omega;

    /* s += alpha·p */
    {
      qpb_complex_double a = {alpha, 0.};
      qpb_spinor_axpy(s, a, p, s);
    }

    /* r -= alpha·D†D·p = alpha·y */
    {
      qpb_complex_double a = {-alpha, 0.};
      qpb_spinor_axpy(r, a, y, r);
    }

    /* new_res_norm = ||r_{k+1}||^2 */
    qpb_spinor_xdotx(&new_res_norm, r);

    /* beta = ||r_{k+1}||^2 / ||r_k||^2 */
    beta = new_res_norm / res_norm;

    /* p_{k+1} = r_{k+1} + beta·p_k */
    {
      qpb_complex_double b_c = {beta, 0.};
      qpb_spinor_axpy(p, b_c, p, r);
    }

    res_norm = new_res_norm;
  }

  if(iters == prec_CG_max_iter)
    return -1;

  return iters;
}


int
qpb_congrad_overlap_kl_pfrac(qpb_spinor_field x, qpb_spinor_field b,
                                        qpb_double CG_epsilon, int CG_max_iter)
{
  qpb_spinor_field p = ov_temp_vecs[9];
  qpb_spinor_field r = ov_temp_vecs[10];
  qpb_spinor_field z = ov_temp_vecs[11];
  qpb_spinor_field y = ov_temp_vecs[12];
  qpb_spinor_field w = ov_temp_vecs[13];
  qpb_spinor_field bprime = ov_temp_vecs[14];
  qpb_spinor_field s = ov_temp_vecs[15];

  int n_reeval = 100;
  int n_echo = 100;
  int iters = 0;

  qpb_double res_norm, true_res_norm, b_norm, bprime_norm;
  qpb_complex_double alpha = {1, 0}, omega = {1, 0};
  qpb_complex_double beta, gamma, new_gamma;

  /* ||b||^2 */
  qpb_spinor_xdotx(&b_norm, b);
  true_res_norm = b_norm;

  /* b' = D†·b */
  qpb_conjugate_overlap_kl_pfrac(bprime, b);
  qpb_spinor_xdotx(&bprime_norm, bprime);

  /* x0 = 0 */
  qpb_spinor_field_set_zero(x);

  /* r0 = b,  z0 = b'  (exact since x0 = 0) */
  qpb_spinor_xeqy(r, b);
  qpb_spinor_xeqy(z, bprime);

  /* Solve M·s0 = z0 for s0 */
  qpb_preconditioner_CG(s, z);

  /* gamma_0 = z0†·s0  (real by HPD of M; take .re explicitly) */
  qpb_spinor_xdoty(&gamma, z, s);
  gamma.im = 0.;

  /* p0 = s0  (preconditioned initial search direction) */
  qpb_spinor_xeqy(p, s);

  qpb_double t = qpb_stop_watch(0);
  for(iters = 1; iters < CG_max_iter; iters++)
  {
    /* Stopping criterion on true residual of original system D·x = b */
    if(true_res_norm / b_norm <= CG_epsilon)
      break;

    /* w = D·p,  y = D†·w  (= A·p) */
    qpb_overlap_kl_pfrac(w, p);
    qpb_conjugate_overlap_kl_pfrac(y, w);

    /* omega = p†·A·p */
    qpb_spinor_xdoty(&omega, p, y);

    /* alpha = gamma / omega */
    alpha = CDEV(gamma, omega);

    /* x += alpha·p */
    qpb_spinor_axpy(x, alpha, p, x);

    /* Update r and z: full recomputation or recursive */
    if(iters % n_reeval == 0)
    {
      qpb_overlap_kl_pfrac(w, x);
      qpb_spinor_xmy(r, b, w);
      qpb_conjugate_overlap_kl_pfrac(y, w);
      qpb_spinor_xmy(z, bprime, y);
    }
    else
    {
      alpha.re = -CDEVR(gamma, omega);
      alpha.im = -CDEVI(gamma, omega);
      qpb_spinor_axpy(r, alpha, w, r);   /* r -= alpha·(D·p)    */
      qpb_spinor_axpy(z, alpha, y, z);   /* z -= alpha·(A·p)    */
    }

    /* Solve M·s = z  (one preconditioner application per outer iteration) */
    qpb_preconditioner_CG(s, z);

    /* new_gamma = z†·s  (real by HPD of M; take .re explicitly) */
    qpb_spinor_xdoty(&new_gamma, z, s);
    res_norm = new_gamma.re;

    /* beta = gamma_{k+1} / gamma_k  (real) */
    beta.re = res_norm / gamma.re;
    beta.im = 0.;

    /* p_{k+1} = s_{k+1} + beta·p_k */
    qpb_spinor_axpy(p, beta, p, s);

    /* Advance gamma */
    gamma.re = res_norm;
    gamma.im = 0.;

    qpb_spinor_xdotx(&true_res_norm, r);
    if(iters % n_echo == 0)
      print(" \t iters = %8d, res = %e\n", iters, true_res_norm / b_norm);
  }
  t = qpb_stop_watch(t);

  /* Final explicit residual check */
  qpb_overlap_kl_pfrac(y, x);
  qpb_spinor_xmy(r, b, y);
  qpb_spinor_xdotx(&true_res_norm, r);

  if(iters == CG_max_iter)
  {
    error(" !\n");
    error(" Preconditioned CG *did not* converge, after %d iterations\n", iters);
    error(" residual = %e, relative = %e, t = %g sec\n", true_res_norm,
                                                          true_res_norm / b_norm, t);
    error(" !\n");
    return -1;
  }

  print(" \tAfter %d iters, preconditioned CG converged, res = %e, relative = %e, "
        "t = %g sec\n",
         iters, true_res_norm, true_res_norm / b_norm, t);

  return iters;
}
