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

static qpb_double rho_plus;
static qpb_double rho_minus;

static qpb_double *numerators;
static qpb_double *shifts;
static qpb_double constant_term;

static qpb_double *even_shifts;
static qpb_double *odd_shifts;

static qpb_double *even_coefficients;
static qpb_double *odd_coefficients;


/* --------------------------- SCALAR FUNCTIONS --------------------------- */ 

/* c_k = tan^2(pi/2 * k/(2n+1)) */
static double ck(int k, int n)
{
    double t = tan((M_PI / 2.0) * ((double)k / (2*n + 1)));
    return t * t;
}

/*
 * Fills the four arrays (all length n).  Caller must free them.
 *
 *  even_shifts[j]  = c_{2(j+1)}       j = 0..n-1   (poles of even product)
 *  even_coeffs[j]  = e_{j+1}
 *
 *  odd_shifts[j]   = c_{2(j+1)-1}     j = 0..n-1   (poles of odd product)
 *  odd_coeffs[j]   = o_{j+1}
 *
 * Even expression:  prod (x^2+c_{2i-1})/(x^2+c_{2i})
 *                 = 1 + sum_j  even_coeffs[j] / (x^2 + even_shifts[j])
 *
 * Odd expression:   prod (x^2+c_{2i})/(x^2+c_{2i-1})
 *                 = 1 + sum_j  odd_coeffs[j]  / (x^2 + odd_shifts[j])
 */
int compute_coefficients(int n,
                         double **even_shifts, double **even_coeffs,
                         double **odd_shifts,  double **odd_coeffs)
{
    *even_shifts = malloc(n * sizeof(double));
    *even_coeffs = malloc(n * sizeof(double));
    *odd_shifts  = malloc(n * sizeof(double));
    *odd_coeffs  = malloc(n * sizeof(double));

    if (!*even_shifts || !*even_coeffs || !*odd_shifts || !*odd_coeffs) {
        free(*even_shifts); free(*even_coeffs);
        free(*odd_shifts);  free(*odd_coeffs);
        return -1;
    }

    /* Pre-cache all c_k values we need (k = 1 .. 2n) */
    double *c = malloc((2*n + 1) * sizeof(double));
    if (!c) return -1;
    for (int k = 1; k <= 2*n; k++)
        c[k] = ck(k, n);

    for (int j = 1; j <= n; j++) {

        /* ---- even: pole at x^2 = -c_{2j} ---- */
        double pole_e = c[2*j];
        (*even_shifts)[j-1] = pole_e;

        double e = c[2*j - 1] - c[2*j];        /* i == j factor */
        for (int i = 1; i <= n; i++) {
            if (i == j) continue;
            e *= (c[2*i - 1] - pole_e) / (c[2*i] - pole_e);
        }
        (*even_coeffs)[j-1] = e;

        /* ---- odd: pole at x^2 = -c_{2j-1} ---- */
        double pole_o = c[2*j - 1];
        (*odd_shifts)[j-1] = pole_o;

        double o = c[2*j] - c[2*j - 1];        /* i == j factor */
        for (int i = 1; i <= n; i++) {
            if (i == j) continue;
            o *= (c[2*i] - pole_o) / (c[2*i - 1] - pole_o);
        }
        (*odd_coeffs)[j-1] = o;
    }

    free(c);
    return 0;
}


/* --------------------------- MATRIX FUNCTIONS --------------------------- */ 

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

    rho_plus = ov_params.rho + 0.5*ov_params.mass;
    rho_minus = ov_params.rho - 0.5*ov_params.mass;

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
      // odd_shifts[i] = pow(tan(trig_arg), 2);
      // even_shifts[i] = 1.0/pow(tan(trig_arg), 2);
      // print("odd_shift[%d] = %.25f, even_shift[%d] = %.25f\n", i, odd_shifts[i], i, even_shifts[i]);
    }

    odd_shifts = qpb_alloc(sizeof(qpb_double)*KL_diagonal_order);
    even_shifts = qpb_alloc(sizeof(qpb_double)*KL_diagonal_order);

    even_coefficients = qpb_alloc(sizeof(qpb_double)*KL_diagonal_order);
    odd_coefficients = qpb_alloc(sizeof(qpb_double)*KL_diagonal_order);

    compute_coefficients(KL_diagonal_order, &even_shifts, &even_coefficients,\
                                               &odd_shifts, &odd_coefficients);

    for(int i=0; i<KL_diagonal_order; i++)
    {
      print("odd_shift[%d] = %.25f, even_shift[%d] = %.25f\n",\
                                           i, odd_shifts[i], i, even_shifts[i]);
      print("odd_coefficients[%d] = %.25f, even_coefficients[%d] = %.25f\n",\
                               i, odd_coefficients[i], i, even_coefficients[i]);
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
  free(odd_shifts);
  free(even_shifts);
  
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
X_op(qpb_spinor_field y, qpb_spinor_field x)
{
  /* Implements X ≡ γ5(a*D - rho*I) */

  void *dslash_args[4];

  dslash_args[0] = ov_params.gauge_ptr;
  dslash_args[1] = &ov_params.m_bare;
  dslash_args[2] = &ov_params.clover;
  dslash_args[3] = &ov_params.c_sw;

  ov_params.g5_dslash_op(y, x, dslash_args);

  return;
}


INLINE void
X_op_shifted(qpb_spinor_field y, qpb_spinor_field x, qpb_double shift)
{
  /* Implements X ≡ γ5(a*D - rho*I) */

  void *dslash_args[4];

  dslash_args[0] = ov_params.gauge_ptr;
  dslash_args[1] = &ov_params.m_bare;
  dslash_args[2] = &ov_params.clover;
  dslash_args[3] = &ov_params.c_sw;

  qpb_spinor_field z = ov_temp_vecs[15];

  ov_params.g5_dslash_op(z, x, dslash_args);

  qpb_spinor_axpy(y, (qpb_complex) {shift, 0.0}, x, z);

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

  qpb_complex a = {rho_plus, 0.};
  qpb_complex b = {rho_minus, 0.};

  qpb_gamma5_sign_function_of_X_pfrac(z, x);

  qpb_spinor_axpby(y, a, x, b, z);

  return;
}


void
qpb_odd_partial_fraction_decomposition(qpb_spinor_field y, qpb_spinor_field x)
{
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
    KL_diagonal_order, odd_shifts, ov_params.c_sw, MS_solver_precision, \
    MS_maximum_solver_iterations);

  // Initialize sum with x
  qpb_spinor_xeqy(sum, x);
  // And then add the rest of the partial fraction terms
  for(int sigma=0; sigma<KL_diagonal_order; sigma++)
    qpb_spinor_axpy(sum, (qpb_complex) {odd_coefficients[sigma], 0.}, yMS[sigma], sum);

  qpb_spinor_xeqy(y, sum);

  return;
}


void
qpb_even_partial_fraction_decomposition(qpb_spinor_field y, qpb_spinor_field x)
{
  qpb_spinor_field sum = ov_temp_vecs[3];

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
    KL_diagonal_order, even_shifts, ov_params.c_sw, MS_solver_precision, \
    MS_maximum_solver_iterations);

  // Initialize sum with x
  qpb_spinor_xeqy(sum, x);
  // And then add the rest of the partial fraction terms
  for(int sigma=0; sigma<KL_diagonal_order; sigma++)
    qpb_spinor_axpy(sum, (qpb_complex) {even_coefficients[sigma], 0.}, yMS[sigma], sum);

  qpb_spinor_xeqy(y, sum);

  return;
}


void
qpb_preconditioner(qpb_spinor_field y, qpb_spinor_field x)
{
  /* Implements:
        Dov,m(x) = 
  */
  
  qpb_spinor_field z = ov_temp_vecs[4];
  qpb_spinor_field w = ov_temp_vecs[5];

  qpb_complex a = {1.0/rho_plus, 0.};
  qpb_complex b = {-1.0/(rho_minus*constant_term), 0.};

  // First term
  X_op(y, x);
  X_op_shifted(z, y, even_shifts[0]);

  // Second term
  qpb_spinor_gamma5(y, x);
  // qpb_even_partial_fraction_decomposition(w, y);
  X_op_shifted(w, y, odd_shifts[0]);

  qpb_spinor_axpby(y, a, z, b, w);

  return;
}


void
qpb_preconditioned_overlap_kl_pfrac(qpb_spinor_field y, qpb_spinor_field x)
{
  /* Implements:
        Dov,m(x) = 
  */
  
  qpb_spinor_field z = ov_temp_vecs[6];
  qpb_spinor_field w = ov_temp_vecs[7];

  qpb_complex a = {rho_minus*constant_term/rho_plus, 0.};
  qpb_complex b = {-rho_plus/(rho_minus*constant_term), 0.};

  // First term
  qpb_odd_partial_fraction_decomposition(y, x);
  D_op(z, y);
  X_op(y, z);
  // Shifted X operator with the second shift of the odd partial fraction expansion
  X_op_shifted(z, y, even_shifts[0]);


  // Second term
  qpb_spinor_gamma5(y, x);
  // qpb_even_partial_fraction_decomposition(w, y);
  X_op_shifted(w, y, odd_shifts[0]);

  qpb_spinor_axpby(y, a, z, b, w);

  return;
}


void
qpb_conjugate_preconditioned_overlap_kl_pfrac(qpb_spinor_field y, qpb_spinor_field x)
{
  /* Implements:
        Dov,m(x) = 
  */
  
  qpb_spinor_field z = ov_temp_vecs[8];
  qpb_spinor_field w = ov_temp_vecs[9];

  qpb_complex a = {rho_minus*constant_term/rho_plus, 0.};
  qpb_complex b = {-rho_plus/(rho_minus*constant_term), 0.};

  // First term
  X_op_shifted(y, x, even_shifts[0]);
  D_op(z, y);
  qpb_odd_partial_fraction_decomposition(y, z);
  X_op(z, y);

  // Second term
  // qpb_even_partial_fraction_decomposition(y, x);
  X_op_shifted(y, x, odd_shifts[0]);
  qpb_spinor_gamma5(w, y);

  qpb_spinor_axpby(y, a, z, b, w);

  return;
}


int
qpb_congrad_overlap_kl_pfrac(qpb_spinor_field x, qpb_spinor_field b, \
                                        qpb_double CG_epsilon, int CG_max_iter)
{
  qpb_spinor_field p = ov_temp_vecs[10];
  qpb_spinor_field r = ov_temp_vecs[11];
  qpb_spinor_field y = ov_temp_vecs[12];
  qpb_spinor_field w = ov_temp_vecs[13];
  qpb_spinor_field bprime = ov_temp_vecs[14];

  int n_reeval = 100;
  int n_echo = 100;
  int iters = 0;
  /**
   * Indicates whether the final inversion test should be performed.
   * Set to 1 (true) to enable the test, or 0 (false) to disable it.
   */
  int final_inversion_test = 1;
  
  qpb_double res_norm, b_norm;
  qpb_complex_double alpha = {1, 0}, omega = {1, 0};
  qpb_complex_double beta, gamma;

  qpb_preconditioner(w, b);
  qpb_conjugate_preconditioned_overlap_kl_pfrac(bprime, w);

  qpb_spinor_xdotx(&b_norm, bprime);

  qpb_spinor_field_set_zero(x);

  /* r0 = bprime - A(x) */
  // qpb_preconditioned_overlap_kl_pfrac(w, x);
  // qpb_conjugate_preconditioned_overlap_kl_pfrac(p, w);
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
    qpb_preconditioned_overlap_kl_pfrac(w, p);
    qpb_conjugate_preconditioned_overlap_kl_pfrac(y, w);

    /* omega = dot(p, A(p)) */
    qpb_spinor_xdoty(&omega, p, y);

    /* alpha = dot(r, r)/omega */
    alpha = CDEV(gamma, omega);

    /* x <- x + alpha*p */
    qpb_spinor_axpy(x, alpha, p, x);

    if(iters % n_reeval == 0) 
    {
      qpb_preconditioned_overlap_kl_pfrac(w, x);
      qpb_conjugate_preconditioned_overlap_kl_pfrac(y, w);
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

  qpb_preconditioned_overlap_kl_pfrac(w, x);
  qpb_conjugate_preconditioned_overlap_kl_pfrac(y, w);
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
  
  if (final_inversion_test)
  {
    qpb_overlap_kl_pfrac(y, x);
    qpb_spinor_xmy(r, b, y);
    qpb_spinor_xdotx(&res_norm, r);
    qpb_spinor_xdotx(&b_norm, b);
    print(" \tFinal relative residual = %e\n", res_norm / b_norm);
  }

  return iters;
}
