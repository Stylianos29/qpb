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


#define OVERLAP_NUMB_TEMP_VECS 19
#define MSCG_NUMB_TEMP_VECS 20


static qpb_spinor_field ov_temp_vecs[OVERLAP_NUMB_TEMP_VECS];
static qpb_spinor_field mscg_temp_vecs[MSCG_NUMB_TEMP_VECS];

static qpb_overlap_params ov_params;

static int KL_diagonal_order;
static qpb_double MS_solver_precision;
static int MS_maximum_solver_iterations;

static qpb_double preconditioner_mass;
static qpb_double prec_CG_epsilon;
static int prec_CG_max_iter;

static qpb_double *numerators;
static qpb_double *shifts;
static qpb_double constant_term;


void
qpb_overlap_kl_pfrac_init(void * gauge, qpb_clover_term clover, \
          enum qpb_kl_classes kl_class, int kl_iters, qpb_double rho, \
          qpb_double c_sw, qpb_double mass, qpb_double scaling_factor, \
          qpb_double ms_epsilon, int ms_max_iters,
          qpb_double prec_mass, qpb_double prec_epsilon, int prec_max_iters)
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
    MS_maximum_solver_iterations = ms_max_iters;

    prec_CG_epsilon = prec_epsilon;
    prec_CG_max_iter = prec_max_iters;
    preconditioner_mass = prec_mass;

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
X_op(qpb_spinor_field y, qpb_spinor_field x)
{
  /* Implements: X ≡ γ5(a*D - ρ) . */

  void *dslash_args[4];

  dslash_args[0] = ov_params.gauge_ptr;
  dslash_args[1] = &ov_params.m_bare;
  dslash_args[2] = &ov_params.clover;
  dslash_args[3] = &ov_params.c_sw;

  ov_params.g5_dslash_op(y, x, dslash_args);

  return;
}


INLINE void
X2_shifted_op(qpb_spinor_field y, qpb_spinor_field x, qpb_double shift)
{
  /* Implements: X^2 + shift . */

  qpb_spinor_field z = ov_temp_vecs[0];

  void *dslash_args[4];

  dslash_args[0] = ov_params.gauge_ptr;
  dslash_args[1] = &ov_params.m_bare;
  dslash_args[2] = &ov_params.clover;
  dslash_args[3] = &ov_params.c_sw;

  ov_params.g5_dslash_op(y, x, dslash_args);
  ov_params.g5_dslash_op(z, y, dslash_args);

  qpb_spinor_axpy(y, (qpb_complex) {shift, 0.0}, x, z);

  return;
}


// INLINE void
// M_op(qpb_spinor_field y, qpb_spinor_field x)
// {
//   /* Implements: M = 3*c*(X^2+1/3) + r*γ5.X(X^2+3) ,
//       with c = ρ + a*m/2, the preconditioner centerand, and
//       r = ρ - a*m/2, the preconditioner radius. */

//   qpb_spinor_field z = ov_temp_vecs[1];
//   qpb_spinor_field w = ov_temp_vecs[2];

//   // qpb_double overlap_mass = ov_params.mass;
//   qpb_double rho = ov_params.rho;

//   qpb_complex three_times_c = {3*(rho + 0.5*preconditioner_mass), 0.};
//   qpb_complex r = {rho - 0.5*preconditioner_mass, 0.};

//   // Part 1
//   X2_shifted_op(z, x, 1.0/3.0);

//   // Part 2
//   X2_shifted_op(y, x, 3.0);
//   D_op(w, y); // γ5.X

//   qpb_spinor_axpby(y, three_times_c, z, r, w);
  
//   return;
// }


// INLINE void
// M_conj_op(qpb_spinor_field y, qpb_spinor_field x)
// {
//   /* Implements: M = 3*c*(X^2+1/3) + r*X(X^2+3)γ5 ,
//       with c = ρ + a*m/2, the preconditioner centerand, and
//       r = ρ - a*m/2, the preconditioner radius. */

//   qpb_spinor_field z = ov_temp_vecs[3];
//   qpb_spinor_field w = ov_temp_vecs[4];

//   // qpb_double overlap_mass = ov_params.mass;
//   qpb_double rho = ov_params.rho;

//   qpb_complex three_times_c = {3*(rho + 0.5*preconditioner_mass), 0.};
//   qpb_complex r = {rho - 0.5*preconditioner_mass, 0.};

//   // Part 1
//   X2_shifted_op(z, x, 1.0/3.0);
  
//   // Part 2
//   qpb_spinor_gamma5(w, x);
//   X2_shifted_op(y, w, 3.0);
//   X_op(w, y);

//   qpb_spinor_axpby(y, three_times_c, z, r, w);
  
//   return;
// }


INLINE void
M_nonsq_op(qpb_spinor_field y, qpb_spinor_field x)
{
  /* Implements the non-squared multiply-up operator for n=1:
        M'' = 3*C*(X^2+1/3)*γ5 + R*X*(X^2+3) ,
     obtained by left-multiplying D_{ov,n=1} by Q_{11}(X^2)*γ5 = 3(X^2+1/3)*γ5.
     This avoids the squaring (M†M) needed in the CG path, since BiCGStab
     handles non-Hermitian systems directly. */

  qpb_spinor_field z = ov_temp_vecs[1];
  qpb_spinor_field w = ov_temp_vecs[2];

  qpb_double rho = ov_params.rho;

  qpb_complex three_times_c = {3*(rho + 0.5*preconditioner_mass), 0.};
  qpb_complex r = {rho - 0.5*preconditioner_mass, 0.};

  /* Part 1:  3C*(X^2+1/3)*γ5*x
     First γ5, then (X^2+1/3): order matters since X^2 does not commute with γ5. */
  qpb_spinor_gamma5(w, x);           /* w  = γ5·x                      */
  X2_shifted_op(z, w, 1.0/3.0);      /* z  = (X^2 + 1/3)·γ5·x          */

  /* Part 2:  R*X*(X^2+3)*x */
  X2_shifted_op(y, x, 3.0);          /* y  = (X^2 + 3)·x               */
  X_op(w, y);                         /* w  = X·(X^2 + 3)·x             */

  /* Combine: y = 3C·z + R·w */
  qpb_spinor_axpby(y, three_times_c, z, r, w);

  return;
}


void
qpb_gamma5_sign_function_of_X_pfrac(qpb_spinor_field y, qpb_spinor_field x)
{
  /* Implements: γ5(sign(X(x))) = γ5(X(c_0 + Sum_{i=1}^{n} c_i/(X^2+σ_i) )),
      with X(x) = γ5(D(x) - ρ*x) . */

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
  
  qpb_spinor_field z = ov_temp_vecs[4];

  qpb_double overlap_mass = ov_params.mass;
  qpb_double rho = ov_params.rho;

  qpb_complex rho_plus = {rho + 0.5*overlap_mass, 0.};
  qpb_complex rho_minus = {rho - 0.5*overlap_mass, 0.};

  qpb_gamma5_sign_function_of_X_pfrac(z, x);

  qpb_spinor_axpby(y, rho_plus, x, rho_minus, z);

  return;
}


// void
// qpb_conjugate_overlap_kl_pfrac(qpb_spinor_field y, qpb_spinor_field x)
// {
//   qpb_spinor_field z = ov_temp_vecs[9];

//   qpb_spinor_gamma5(y, x);
//   qpb_overlap_kl_pfrac(z, y);
//   qpb_spinor_gamma5(y, z);

//   return;
// }


// int
// qpb_preconditioner_CG(qpb_spinor_field x, qpb_spinor_field b)
// {
//   qpb_spinor_field r = ov_temp_vecs[10];
//   qpb_spinor_field p = ov_temp_vecs[11];
//   qpb_spinor_field w = ov_temp_vecs[12];
//   qpb_spinor_field y = ov_temp_vecs[13];
//   qpb_spinor_field bprime = ov_temp_vecs[14];
//   qpb_spinor_field s = ov_temp_vecs[15];


//   int iters = 0;

//   /* All scalars are real -- D†D is HPD */
//   qpb_double b_norm, b_prime_norm, res_norm, new_res_norm, omega, alpha, beta;

//   /* ||b||^2 */
//   qpb_spinor_xdotx(&b_norm, b);

//   /* b_prime = 3(X^2+1/3) b */ 
//   X2_shifted_op(w, b, 1.0/3.0);
//   qpb_spinor_ax(bprime, (qpb_complex) {3.0, 0.0}, w);

//   /* s0 = 0,  r0 = b_prime */
//   qpb_spinor_field_set_zero(s);
//   qpb_spinor_xeqy(r, bprime);
//   qpb_spinor_xdotx(&b_prime_norm, bprime);

//   /* gamma_0 = ||r0||^2 */
//   qpb_spinor_xdotx(&res_norm, r);

//   /* p0 = r0 */
//   qpb_spinor_xeqy(p, r);

//   for(iters = 1; iters < prec_CG_max_iter; iters++)
//   {
//     /* Stopping on relative residual of the HPD system:
//        ||r||^2 / ||b||^2 <= prec_CG_epsilon^2   (squaring avoids a sqrt) */
//     if(res_norm / b_prime_norm <= prec_CG_epsilon * prec_CG_epsilon)
//       break;

//     /* Apply D†D to p: w = D·p,  y = D†·w = D†D·p */
//     M_op(w, p);
//     M_conj_op(y, w);

//     /* omega = p†·D†D·p = ||D·p||^2 = ||w||^2  (real, positive) */
//     qpb_spinor_xdotx(&omega, w);

//     /* alpha = ||r||^2 / (p†·D†D·p)  (real, positive) */
//     alpha = res_norm / omega;

//     /* s += alpha·p */
//     {
//       qpb_complex a = {alpha, 0.};
//       qpb_spinor_axpy(s, a, p, s);
//     }

//     /* r -= alpha·D†D·p = alpha·y */
//     {
//       qpb_complex a = {-alpha, 0.};
//       qpb_spinor_axpy(r, a, y, r);
//     }

//     /* new_res_norm = ||r_{k+1}||^2 */
//     qpb_spinor_xdotx(&new_res_norm, r);

//     /* beta = ||r_{k+1}||^2 / ||r_k||^2 */
//     beta = new_res_norm / res_norm;

//     /* p_{k+1} = r_{k+1} + beta·p_k */
//     {
//       qpb_complex b_c = {beta, 0.};
//       qpb_spinor_axpy(p, b_c, p, r);
//     }

//     res_norm = new_res_norm;
//   }

//   X2_shifted_op(y, s, 1.0/3.0);
//   qpb_spinor_ax(x, (qpb_complex) {3.0, 0.0}, y);

//   if(iters == prec_CG_max_iter)
//     return -1;

//   return iters;
// }


int
qpb_preconditioner_bicgstab(qpb_spinor_field x, qpb_spinor_field b)
{
  /* Solves the non-squared multiply-up system  M''·x = b'  via BiCGStab,
     where  b' = 3(X^2+1/3)·γ5·b  and  M'' = 3C(X^2+1/3)γ5 + R·X(X^2+3).
     The solution x satisfies  D_{ov,n=1}·x = b  directly (no post-processing).

     This replaces qpb_preconditioner_CG for the BiCGStab path: since we
     no longer need D†D, there is no squaring — one fewer operator application
     per inner iteration. */

  qpb_spinor_field r0_hat = ov_temp_vecs[5];
  qpb_spinor_field r      = ov_temp_vecs[6];
  qpb_spinor_field p      = ov_temp_vecs[7];
  qpb_spinor_field u      = ov_temp_vecs[8];
  qpb_spinor_field v      = ov_temp_vecs[9];
  qpb_spinor_field bprime = ov_temp_vecs[10];
  qpb_spinor_field tmp    = ov_temp_vecs[11]; /* scratch for b' computation */

  int iters = 0;
  const int n_reeval = 10;

  qpb_double res_norm, bprime_norm;
  qpb_complex_double alpha = {1, 0}, omega = {1, 0};
  qpb_complex_double beta, gamma, rho, zeta;

  /* --- Compute modified RHS:  b' = 3(X^2+1/3)·γ5·b --- */
  qpb_spinor_gamma5(r, b);                 /* r    = γ5·b - r as temp here   */
  X2_shifted_op(tmp, r, 1.0/3.0);       /* tmp = (X^2+1/3)·γ5·b      */
  qpb_spinor_ax(bprime, (qpb_complex) {3.0, 0.0}, tmp);
                                              /* bprime = 3(X^2+1/3)·γ5·b     */
  qpb_spinor_xdotx(&bprime_norm, bprime);

  /* --- x0 = 0 --- */
  qpb_spinor_field_set_zero(x);
  qpb_spinor_field_set_zero(p);
  qpb_spinor_field_set_zero(u);
  qpb_spinor_field_set_zero(v);

  /* --- r0 = b' - M''·x0 = b'  (since x0 = 0) --- */
  qpb_spinor_xeqy(r, bprime);
  qpb_spinor_xeqy(r0_hat, r);        /* shadow residual (fixed) */

  qpb_spinor_xdotx(&gamma.re, r);
  gamma.im = 0;
  res_norm = gamma.re;
  rho = gamma;

  for(iters = 1; iters < prec_CG_max_iter; iters++)
  {
    print(" \tAt %d iters, preconditioner solver has, " \
                  "relative res = %e\n", iters, res_norm / bprime_norm);

    if(res_norm / bprime_norm <= prec_CG_epsilon * prec_CG_epsilon)
      break;

    /* gamma = (r̂0, r) */
    qpb_spinor_xdoty(&gamma, r0_hat, r);

    /* beta = (gamma/rho)·(alpha/omega) */
    beta = CMUL(CDEV(gamma, rho), CDEV(alpha, omega));

    /* p = r + beta·(p - omega·u) */
    omega = CNEGATE(omega);
    qpb_spinor_axpy(p, omega, u, p);      /* p -= omega·u                */
    qpb_spinor_axpy(p, beta, p, r);       /* p  = beta·p + r             */

    /* u = M''·p */
    M_nonsq_op(u, p);

    /* rho = gamma,  alpha = rho / (r̂0, u) */
    qpb_spinor_xdoty(&beta, r0_hat, u);
    rho = gamma;
    alpha = CDEV(rho, beta);

    /* r -= alpha·u */
    alpha = CNEGATE(alpha);
    qpb_spinor_axpy(r, alpha, u, r);

    /* v = M''·r  (= M''·s, the 't' vector) */
    M_nonsq_op(v, r);

    /* omega = (v, r) / (v, v) */
    qpb_spinor_xdoty(&zeta, v, r);
    qpb_spinor_xdotx(&beta.re, v);
    beta.im = 0;
    omega = CDEV(zeta, beta);

    /* x += alpha·p + omega·r  (recall alpha is currently negated) */
    alpha = CNEGATE(alpha);
    qpb_spinor_axpy(x, omega, r, x);     /* x += omega·s               */
    qpb_spinor_axpy(x, alpha, p, x);     /* x += alpha·p               */

    /* Residual update */
    if(iters % n_reeval == 0)
    {
      M_nonsq_op(r, x);
      qpb_spinor_xmy(r, bprime, r);      /* r = b' - M''·x             */
    }
    else
    {
      omega = CNEGATE(omega);
      qpb_spinor_axpy(r, omega, v, r);   /* r -= omega·v               */
      omega = CNEGATE(omega);
    }

    qpb_spinor_xdotx(&res_norm, r);
  }

  if(iters == prec_CG_max_iter)
    return -1;


  return iters;
}


// int
// qpb_congrad_overlap_kl_pfrac(qpb_spinor_field x, qpb_spinor_field b,
//                                         qpb_double CG_epsilon, int CG_max_iter)
// {
//   qpb_spinor_field p = ov_temp_vecs[23];
//   qpb_spinor_field r = ov_temp_vecs[24];
//   qpb_spinor_field z = ov_temp_vecs[25];
//   qpb_spinor_field y = ov_temp_vecs[26];
//   qpb_spinor_field w = ov_temp_vecs[27];
//   qpb_spinor_field bprime = ov_temp_vecs[28];
//   qpb_spinor_field s = ov_temp_vecs[29];

//   int n_reeval = 100;
//   int n_echo = 100;
//   int iters = 0;

//   qpb_double res_norm, true_res_norm, b_norm, bprime_norm;
//   qpb_complex alpha = {1, 0}, omega = {1, 0};
//   qpb_complex beta, gamma, new_gamma;

//   /* ||b||^2 */
//   qpb_spinor_xdotx(&b_norm, b);
//   true_res_norm = b_norm;

//   /* b' = D†·b */
//   qpb_conjugate_overlap_kl_pfrac(bprime, b);
//   qpb_spinor_xdotx(&bprime_norm, bprime);

//   /* x0 = 0 */
//   qpb_spinor_field_set_zero(x);

//   /* r0 = b,  z0 = b'  (exact since x0 = 0) */
//   qpb_spinor_xeqy(r, b);
//   qpb_spinor_xeqy(z, bprime);

//   /* Solve M·s0 = z0 for s0 (or copy if no preconditioning) */
//   if(prec_CG_max_iter == 0)
//     qpb_spinor_xeqy(s, z);
//   else
//     qpb_preconditioner_CG(s, z);

//   /* gamma_0 = z0†·s0  (real by HPD of M; take .re explicitly) */
//   qpb_spinor_xdoty(&gamma, z, s);
//   gamma.im = 0.;

//   /* p0 = s0  (preconditioned initial search direction) */
//   qpb_spinor_xeqy(p, s);

//   qpb_double t = qpb_stop_watch(0);
//   for(iters = 1; iters < CG_max_iter; iters++)
//   {
//     /* Stopping criterion on true residual of original system D·x = b */
//     if(true_res_norm / b_norm <= CG_epsilon)
//       break;

//     /* w = D·p,  y = D†·w  (= A·p) */
//     qpb_overlap_kl_pfrac(w, p);
//     qpb_conjugate_overlap_kl_pfrac(y, w);

//     /* omega = p†·A·p */
//     qpb_spinor_xdoty(&omega, p, y);

//     /* alpha = gamma / omega */
//     alpha = CDEV(gamma, omega);

//     /* x += alpha·p */
//     qpb_spinor_axpy(x, alpha, p, x);

//     /* Update r and z: full recomputation or recursive */
//     if(iters % n_reeval == 0)
//     {
//       qpb_overlap_kl_pfrac(w, x);
//       qpb_spinor_xmy(r, b, w);
//       qpb_conjugate_overlap_kl_pfrac(y, w);
//       qpb_spinor_xmy(z, bprime, y);
//     }
//     else
//     {
//       alpha.re = -CDEVR(gamma, omega);
//       alpha.im = -CDEVI(gamma, omega);
//       qpb_spinor_axpy(r, alpha, w, r);   /* r -= alpha·(D·p)    */
//       qpb_spinor_axpy(z, alpha, y, z);   /* z -= alpha·(A·p)    */
//     }

//     /* Solve M·s = z for s (or copy if no preconditioning) */
//     if(prec_CG_max_iter == 0)
//       qpb_spinor_xeqy(s, z);
//     else
//       qpb_preconditioner_CG(s, z);

//     /* new_gamma = z†·s  (real by HPD of M; take .re explicitly) */
//     qpb_spinor_xdoty(&new_gamma, z, s);
//     res_norm = new_gamma.re;

//     /* beta = gamma_{k+1} / gamma_k  (real) */
//     beta.re = res_norm / gamma.re;
//     beta.im = 0.;

//     /* p_{k+1} = s_{k+1} + beta·p_k */
//     qpb_spinor_axpy(p, beta, p, s);

//     /* Advance gamma */
//     gamma.re = res_norm;
//     gamma.im = 0.;

//     qpb_spinor_xdotx(&true_res_norm, r);
//     if(iters % n_echo == 0)
//       print(" \t iters = %8d, res = %e\n", iters, true_res_norm / b_norm);
//   }
//   t = qpb_stop_watch(t);

//   /* Final explicit residual check */
//   qpb_overlap_kl_pfrac(y, x);
//   qpb_spinor_xmy(r, b, y);
//   qpb_spinor_xdotx(&true_res_norm, r);

//   if(iters == CG_max_iter)
//   {
//     error(" !\n");
//     error(" Preconditioned CG *did not* converge, after %d iterations\n", iters);
//     error(" residual = %e, relative = %e, t = %g sec\n", true_res_norm,
//                                                           true_res_norm / b_norm, t);
//     error(" !\n");
//     return -1;
//   }

//   print(" \tAfter %d iters, preconditioned CG converged, res = %e, relative = %e, "
//         "t = %g sec\n",
//          iters, true_res_norm, true_res_norm / b_norm, t);

//   return iters;
// }


int
qpb_bicgstab_overlap_kl_pfrac(qpb_spinor_field x, qpb_spinor_field b,
                               qpb_double epsilon, int max_iter)
{
  /* Preconditioned BiCGStab for the overlap operator.
     Solves  D_ov · x = b  directly (no normal equations, no squaring).

     Preconditioner:  K ≈ D_{ov,n=1}  (overlap at diagonal KL order n=1),
     applied via the non-squared multiply-up trick:
       K·y = r  ⟹  M''·y = 3(X^2+1/3)γ5·r ,
     solved by the inner BiCGStab (qpb_preconditioner_bicgstab).

     If prec_CG_max_iter == 0, no preconditioning is used (K = I). */

  qpb_spinor_field r0_hat = ov_temp_vecs[12];
  qpb_spinor_field r      = ov_temp_vecs[13];
  qpb_spinor_field p      = ov_temp_vecs[14];
  qpb_spinor_field u      = ov_temp_vecs[15];   /* A·y  or  A·p (unprec) */
  qpb_spinor_field v      = ov_temp_vecs[16];   /* A·z  or  A·s (unprec) */
  qpb_spinor_field y      = ov_temp_vecs[17];   /* K^{-1}·p              */
  qpb_spinor_field z_vec  = ov_temp_vecs[18];   /* K^{-1}·s              */

  int n_reeval = 100;
  int n_echo = 100;
  int iters = 0;

  qpb_double res_norm, b_norm;
  qpb_complex_double alpha = {1, 0}, omega = {1, 0};
  qpb_complex_double beta, gamma, rho, zeta;

  /* ||b||^2 */
  qpb_spinor_xdotx(&b_norm, b);

  /* x0 = 0 */
  qpb_spinor_field_set_zero(x);

  /* r0 = b - D_ov·x0 = b */
  qpb_spinor_xeqy(r, b);
  qpb_spinor_xeqy(r0_hat, r);

  qpb_spinor_xdotx(&gamma.re, r);
  gamma.im = 0;
  res_norm = gamma.re;
  rho = gamma;

  qpb_double t = qpb_stop_watch(0);
  for(iters = 1; iters < max_iter; iters++)
  {
    if(res_norm / b_norm <= epsilon)
      break;

    /* gamma = (r̂0, r) */
    qpb_spinor_xdoty(&gamma, r0_hat, r);

    /* beta = (gamma/rho)·(alpha/omega) */
    beta = CMUL(CDEV(gamma, rho), CDEV(alpha, omega));

    /* p = r + beta·(p - omega·u) */
    omega = CNEGATE(omega);
    qpb_spinor_axpy(p, omega, u, p);      /* p -= omega·u                */
    qpb_spinor_axpy(p, beta, p, r);       /* p  = beta·p + r             */

    /* --- Preconditioner:  y = K^{-1}·p --- */
    if(prec_CG_max_iter == 0)
      qpb_spinor_xeqy(y, p);             /* no preconditioning           */
    else
      qpb_preconditioner_bicgstab(y, p);

    /* u = D_ov · y */
    qpb_overlap_kl_pfrac(u, y);

    /* rho = gamma,  alpha = rho / (r̂0, u) */
    qpb_spinor_xdoty(&beta, r0_hat, u);
    rho = gamma;
    alpha = CDEV(rho, beta);

    /* r -= alpha·u  (r now holds 's') */
    alpha = CNEGATE(alpha);
    qpb_spinor_axpy(r, alpha, u, r);

    /* --- Preconditioner:  z = K^{-1}·s  (s is stored in r) --- */
    if(prec_CG_max_iter == 0)
      qpb_spinor_xeqy(z_vec, r);
    else
      qpb_preconditioner_bicgstab(z_vec, r);

    /* v = D_ov · z */
    qpb_overlap_kl_pfrac(v, z_vec);

    /* omega = (v, s) / (v, v)  with s stored in r */
    qpb_spinor_xdoty(&zeta, v, r);
    qpb_spinor_xdotx(&beta.re, v);
    beta.im = 0;
    omega = CDEV(zeta, beta);

    /* x += alpha·y + omega·z  (alpha is currently negated → flip back) */
    alpha = CNEGATE(alpha);
    qpb_spinor_axpy(x, alpha, y, x);     /* x += alpha·y                */
    qpb_spinor_axpy(x, omega, z_vec, x); /* x += omega·z                */

    /* Residual update */
    if(iters % n_reeval == 0)
    {
      qpb_overlap_kl_pfrac(u, x);
      qpb_spinor_xmy(r, b, u);           /* r = b - D_ov·x              */
    }
    else
    {
      omega = CNEGATE(omega);
      qpb_spinor_axpy(r, omega, v, r);   /* r -= omega·v                */
      omega = CNEGATE(omega);
    }

    qpb_spinor_xdotx(&res_norm, r);
    if(iters % n_echo == 0)
      print(" \t iters = %8d, res = %e\n", iters, res_norm / b_norm);
  }
  t = qpb_stop_watch(t);

  /* Final explicit residual check */
  qpb_overlap_kl_pfrac(u, x);
  qpb_spinor_xmy(r, b, u);
  qpb_spinor_xdotx(&res_norm, r);

  if(iters == max_iter)
  {
    error(" !\n");
    error(" Preconditioned BiCGStab *did not* converge, after %d iterations\n",
          iters);
    error(" residual = %e, relative = %e, t = %g sec\n", res_norm,
          res_norm / b_norm, t);
    error(" !\n");
    return -1;
  }

  print(" \tAfter %d iters, preconditioned BiCGStab converged, res = %e, "
        "relative = %e, t = %g sec\n",
        iters, res_norm, res_norm / b_norm, t);

  return iters;
}