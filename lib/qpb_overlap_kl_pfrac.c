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


/*============================================================================*/
/*                       SECTION 1 - statics and defines                      */
/*============================================================================*/

/*
  ov_temp_vecs layout (21 vectors total; budget is 25):

    [ 0]    X2_shifted_op scratch                    (nPrec=1 paths)
            M_kernel_op    D·x scratch               (nPrec=0 paths)

    [ 1]    M_kernel_conj_op  γ5·x scratch           (nPrec=0 paths)
            M_op / M_mult_up_op scratch              (nPrec=1 paths)
    [ 2]    M_op / M_mult_up_op scratch              (nPrec=1 paths)
    [ 3]    M_conj_op / M_mult_up_conj_op scratch    (nPrec=1 paths)
    [ 4]    M_conj_op / M_mult_up_conj_op scratch    (nPrec=1 paths)

    [ 5]    qpb_gamma5_sign_function_of_X_pfrac
    [ 6]    qpb_overlap_kl_pfrac
    [ 7]    qpb_conjugate_overlap_kl_pfrac

    [ 8..13]   inner solver (qpb_preconditioner_CG  OR  qpb_preconditioner_CGNE)
                  r, p, w, y, bprime, s|bpp                              (6 vecs)

    [14..20]   outer solver (qpb_congrad_overlap_kl_pfrac  OR
                             qpb_bicgstab_overlap_kl_pfrac)              (7 vecs)

  Indices [8..13] and [14..20] are reused across CG/CGNE inner and CG/BiCGStab
  outer because outer↔inner pairing is strict (outer CG ↔ inner CG, outer
  BiCGStab ↔ inner CGNE) and only one outer solver runs per invocation.

  -  Section 10  -  LEGACY / EXPERIMENTAL: inner BiCGStab for nPrec ∈ {0, 1}
     (guarded by #if 0)
  */

#define OVERLAP_NUMB_TEMP_VECS 21
#define MSCG_NUMB_TEMP_VECS 20


static qpb_spinor_field ov_temp_vecs[OVERLAP_NUMB_TEMP_VECS];
static qpb_spinor_field mscg_temp_vecs[MSCG_NUMB_TEMP_VECS];

static qpb_overlap_params ov_params;

static int KL_diagonal_order;
static qpb_double MS_solver_precision;
static int MS_maximum_solver_iterations;

static int prec_order;
static qpb_double preconditioner_mass;
static qpb_double prec_CG_epsilon;
static int prec_CG_max_iter;

static qpb_double *numerators;
static qpb_double *shifts;
static qpb_double constant_term;


/*
  Preconditioner-operator function pointers, wired in qpb_overlap_kl_pfrac_init()
  based on prec_order:

    Inner standard CG (outer-CGNE path)         → M_inner_op_ptr, M_inner_conj_op_ptr
      nPrec=0: kernel preconditioner            M_kernel_op,  M_kernel_conj_op
      nPrec=1: squared multiply-up operators    M_op,         M_conj_op

    Inner CGNE (outer-BiCGStab path)            → K_inner_op_ptr, K_inner_conj_op_ptr
      nPrec=0: kernel preconditioner            M_kernel_op,    M_kernel_conj_op
      nPrec=1: multiply-up operators            M_mult_up_op,   M_mult_up_conj_op
*/
static void (*M_inner_op_ptr)     (qpb_spinor_field, qpb_spinor_field);
static void (*M_inner_conj_op_ptr)(qpb_spinor_field, qpb_spinor_field);
static void (*K_inner_op_ptr)     (qpb_spinor_field, qpb_spinor_field);
static void (*K_inner_conj_op_ptr)(qpb_spinor_field, qpb_spinor_field);


/* Forward declarations of the operators assigned to the function pointers */
static void M_kernel_op       (qpb_spinor_field y, qpb_spinor_field x);
static void M_kernel_conj_op  (qpb_spinor_field y, qpb_spinor_field x);
static void M_op              (qpb_spinor_field y, qpb_spinor_field x);
static void M_conj_op         (qpb_spinor_field y, qpb_spinor_field x);
static void M_mult_up_op      (qpb_spinor_field y, qpb_spinor_field x);
static void M_mult_up_conj_op (qpb_spinor_field y, qpb_spinor_field x);

/* Forward declaration for the legacy inner BiCGStab preconditioner
   defined in Section 10.  Guarded by the same marker so the symbol
   only exists when the dead-code section is enabled. */
#if 1  /* LEGACY_BICGSTAB_PREC */
  static int qpb_preconditioner_bicgstab(qpb_spinor_field x, qpb_spinor_field b);
#endif /* LEGACY_BICGSTAB_PREC */


/*============================================================================*/
/*                       SECTION 2 - init and finalize                        */
/*============================================================================*/

void
qpb_overlap_kl_pfrac_init(void * gauge, qpb_clover_term clover, \
          enum qpb_kl_classes kl_class, int kl_iters, qpb_double rho, \
          qpb_double c_sw, qpb_double mass, qpb_double scaling_factor, \
          qpb_double ms_epsilon, int ms_max_iters, \
          int prec_order_arg, qpb_double prec_mass, \
                              qpb_double prec_epsilon, int prec_max_iters)
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

    /* Defensive validation: main.c already enforces 0/1, but a single check
       here keeps the library safe to call from any future caller. */
    if(prec_order_arg != 0 && prec_order_arg != 1)
    {
      error("qpb_overlap_kl_pfrac_init: prec_order must be 0 or 1, got %d\n",
            prec_order_arg);
      exit(QPB_PARAMETERS_ERROR);
    }
    prec_order = prec_order_arg;

    /* Wire preconditioner-operator function pointers based on prec_order */
    switch(prec_order)
    {
    case 0:
      /* nPrec=0: both inner solvers use the kernel preconditioner */
      M_inner_op_ptr      = &M_kernel_op;
      M_inner_conj_op_ptr = &M_kernel_conj_op;
      K_inner_op_ptr      = &M_kernel_op;
      K_inner_conj_op_ptr = &M_kernel_conj_op;
      break;
    case 1:
      /* nPrec=1: inner CG uses the squared M, inner CGNE uses the
         multiply-up M''  */
      M_inner_op_ptr      = &M_op;
      M_inner_conj_op_ptr = &M_conj_op;
      K_inner_op_ptr      = &M_mult_up_op;
      K_inner_conj_op_ptr = &M_mult_up_conj_op;
      break;
    }

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


/*============================================================================*/
/*                       SECTION 3 - basic operators                          */
/*============================================================================*/

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


/*============================================================================*/
/*           SECTION 4 - nPrec=0 preconditioner operators (kernel)            */
/*============================================================================*/

static void
M_kernel_op(qpb_spinor_field y, qpb_spinor_field x)
{
  /* Implements M = c·I + r·D ,
     with c = ρ + a*m/2 (preconditioner center)
     and  r = ρ - a*m/2 (preconditioner radius).
     This is the simple kernel-based preconditioner used for nPrec=0. */

  qpb_spinor_field w = ov_temp_vecs[0];

  qpb_double rho = ov_params.rho;

  qpb_complex preconditioner_center = {rho + 0.5*preconditioner_mass, 0.};
  qpb_complex preconditioner_radius = {rho - 0.5*preconditioner_mass, 0.};

  D_op(w, x);
  qpb_spinor_axpby(y, preconditioner_center, x, preconditioner_radius, w);
  
  return;
}


static void
M_kernel_conj_op(qpb_spinor_field y, qpb_spinor_field x)
{
  /* Implements M† = γ5·M·γ5 ,
     valid because of the γ5-hermiticity of D: γ5·D·γ5 = D†. */

  qpb_spinor_field z = ov_temp_vecs[1];

  qpb_spinor_gamma5(y, x);
  M_kernel_op(z, y);
  qpb_spinor_gamma5(y, z);
  
  return;
}


/*============================================================================*/
/*       SECTION 5 - nPrec=1 squared preconditioner operators (CG path)       */
/*============================================================================*/

static void
M_op(qpb_spinor_field y, qpb_spinor_field x)
{
  /* Implements: M = 3*c*(X^2+1/3) + r*γ5.X(X^2+3) ,
      with c = ρ + a*m/2, the preconditioner centerand, and
      r = ρ - a*m/2, the preconditioner radius. */

  qpb_spinor_field z = ov_temp_vecs[1];
  qpb_spinor_field w = ov_temp_vecs[2];

  qpb_double rho = ov_params.rho;

  qpb_complex three_times_c = {3*(rho + 0.5*preconditioner_mass), 0.};
  qpb_complex r = {rho - 0.5*preconditioner_mass, 0.};

  // Part 1
  X2_shifted_op(z, x, 1.0/3.0);

  // Part 2
  X2_shifted_op(y, x, 3.0);
  D_op(w, y); // γ5.X

  qpb_spinor_axpby(y, three_times_c, z, r, w);

  return;
}


static void
M_conj_op(qpb_spinor_field y, qpb_spinor_field x)
{
  /* Implements: M† = 3*c*(X^2+1/3) + r*X(X^2+3)γ5 ,
      with c = ρ + a*m/2, the preconditioner centerand, and
      r = ρ - a*m/2, the preconditioner radius. */

  qpb_spinor_field z = ov_temp_vecs[3];
  qpb_spinor_field w = ov_temp_vecs[4];

  qpb_double rho = ov_params.rho;

  qpb_complex three_times_c = {3*(rho + 0.5*preconditioner_mass), 0.};
  qpb_complex r = {rho - 0.5*preconditioner_mass, 0.};

  // Part 1
  X2_shifted_op(z, x, 1.0/3.0);

  // Part 2
  qpb_spinor_gamma5(w, x);
  X2_shifted_op(y, w, 3.0);
  X_op(w, y);

  qpb_spinor_axpby(y, three_times_c, z, r, w);

  return;
}


/*============================================================================*/
/*    SECTION 6 - nPrec=1 multiply-up preconditioner operators (CGNE path)    */
/*============================================================================*/

static void
M_mult_up_op(qpb_spinor_field y, qpb_spinor_field x)
{
  /* Implements the "multiply-up" operator for n=1:
        M'' = 3*C*(X^2+1/3)*γ5 + R*X*(X^2+3) ,
     obtained by left-multiplying D_{ov,n=1} by 3(X^2+1/3)*γ5.

     This is the operator used by the inner CGNE (outer-BiCGStab path) for
     nPrec=1. Combined with the modified RHS b' = 3(X^2+1/3)·γ5·b , solving
     M''·y = b' yields y satisfying D_{ov,n=1}·y = b directly. */

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
  X_op(w, y);                        /* w  = X·(X^2 + 3)·x             */

  /* Combine: y = 3C·z + R·w */
  qpb_spinor_axpby(y, three_times_c, z, r, w);

  return;
}


static void
M_mult_up_conj_op(qpb_spinor_field y, qpb_spinor_field x)
{
  /* Implements the conjugate of the multiply-up operator:
        (M'')† = 3*C*γ5*(X^2+1/3) + R*X*(X^2+3) .

     Derivation: (M'')† = [3C(X^2+1/3)γ5]† + [R*X(X^2+3)]†
                       = γ5·(X^2+1/3)·3C  +  (X^2+3)·X·R
                       = 3C·γ5·(X^2+1/3)  +  R·X·(X^2+3) ,
     using γ5† = γ5 , (X^2+α)† = X^2+α , X† = X (γ5-hermiticity of D), and
     the commutativity of X and X^2 within polynomial expressions. */

  qpb_spinor_field z = ov_temp_vecs[3];
  qpb_spinor_field w = ov_temp_vecs[4];

  qpb_double rho = ov_params.rho;

  qpb_complex three_times_c = {3*(rho + 0.5*preconditioner_mass), 0.};
  qpb_complex r = {rho - 0.5*preconditioner_mass, 0.};

  /* Part 1:  3C·γ5·(X^2+1/3)·x
     First (X^2+1/3), then γ5: order is the opposite of M_mult_up_op's Part 1. */
  X2_shifted_op(w, x, 1.0/3.0);      /* w = (X^2 + 1/3)·x          */
  qpb_spinor_gamma5(z, w);           /* z = γ5·(X^2 + 1/3)·x       */

  /* Part 2:  R·X·(X^2+3)·x  (identical to M_mult_up_op's Part 2) */
  X2_shifted_op(y, x, 3.0);          /* y = (X^2 + 3)·x            */
  X_op(w, y);                        /* w = X·(X^2 + 3)·x          */

  /* Combine: y = 3C·z + R·w */
  qpb_spinor_axpby(y, three_times_c, z, r, w);

  return;
}


/*============================================================================*/
/*                   SECTION 7 - sign function and overlap                    */
/*============================================================================*/

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


void
qpb_conjugate_overlap_kl_pfrac(qpb_spinor_field y, qpb_spinor_field x)
{
  qpb_spinor_field z = ov_temp_vecs[7];

  qpb_spinor_gamma5(y, x);
  qpb_overlap_kl_pfrac(z, y);
  qpb_spinor_gamma5(y, z);

  return;
}


/*============================================================================*/
/*                      SECTION 8 - inner solvers                             */
/*============================================================================*/

int
qpb_preconditioner_CG(qpb_spinor_field x, qpb_spinor_field b)
{
  /* Inner standard CG solver invoked by the outer CGNE (qpb_congrad_overlap_kl_pfrac).
     Solves the HPD system (M†M)·s = bprime where M is selected by prec_order:

       nPrec=0:  M = M_kernel_op,    bprime = b,                  x = s
       nPrec=1:  M = M_op,           bprime = 3(X^2+1/3)·b,       x = 3(X^2+1/3)·s

       The nPrec=1 algorithm reproduces the previous (validated) implementation
     bit-for-bit; the nPrec=0 path uses bprime = b and skips post-processing,
     because the kernel preconditioner does not need the multiply-up trick. */

  qpb_spinor_field r      = ov_temp_vecs[8];
  qpb_spinor_field p      = ov_temp_vecs[9];
  qpb_spinor_field w      = ov_temp_vecs[10];
  qpb_spinor_field y      = ov_temp_vecs[11];
  qpb_spinor_field bprime = ov_temp_vecs[12];
  qpb_spinor_field s      = ov_temp_vecs[13];

  int iters = 0;

  /* All scalars are real -- M†M is HPD */
  qpb_double b_prime_norm, res_norm, new_res_norm, omega, alpha, beta;

  /* --- bprime preprocessing --- */
  if(prec_order == 1)
  {
    /* bprime = 3(X^2+1/3)·b */
    X2_shifted_op(w, b, 1.0/3.0);
    qpb_spinor_ax(bprime, (qpb_complex) {3.0, 0.0}, w);
  }
  else /* prec_order == 0 */
  {
    /* bprime = b */
    qpb_spinor_xeqy(bprime, b);
  }
  qpb_spinor_xdotx(&b_prime_norm, bprime);

  /* s0 = 0,  r0 = bprime */
  qpb_spinor_field_set_zero(s);
  qpb_spinor_xeqy(r, bprime);

  /* gamma_0 = ||r0||^2 */
  qpb_spinor_xdotx(&res_norm, r);

  /* p0 = r0 */
  qpb_spinor_xeqy(p, r);

  for(iters = 1; iters < prec_CG_max_iter; iters++)
  {
    /* Stopping on relative residual of the HPD system:
       ||r||^2 / ||bprime||^2 <= prec_CG_epsilon^2   (squaring avoids a sqrt) */
    if(res_norm / b_prime_norm <= prec_CG_epsilon * prec_CG_epsilon)
      break;

    /* Apply M†M to p: w = M·p,  y = M†·w = M†M·p */
    M_inner_op_ptr(w, p);
    M_inner_conj_op_ptr(y, w);

    /* omega = p†·M†M·p = ||M·p||^2 = ||w||^2  (real, positive) */
    qpb_spinor_xdotx(&omega, w);

    /* alpha = ||r||^2 / (p†·M†M·p)  (real, positive) */
    alpha = res_norm / omega;

    /* s += alpha·p */
    qpb_spinor_axpy(s, (qpb_complex) {alpha, 0.}, p, s);

    /* r -= alpha·M†M·p = alpha·y */
    qpb_spinor_axpy(r, (qpb_complex) {-alpha, 0.}, y, r);

    /* new_res_norm = ||r_{k+1}||^2 */
    qpb_spinor_xdotx(&new_res_norm, r);

    /* beta = ||r_{k+1}||^2 / ||r_k||^2 */
    beta = new_res_norm / res_norm;

    /* p_{k+1} = r_{k+1} + beta·p_k */
    qpb_spinor_axpy(p, (qpb_complex) {beta, 0.}, p, r);

    res_norm = new_res_norm;
  }
  
  /* --- postprocessing: recover x from s --- */
  if(prec_order == 1)
  {
    /* x = 3(X^2+1/3)·s */
    X2_shifted_op(y, s, 1.0/3.0);
    qpb_spinor_ax(x, (qpb_complex) {3.0, 0.0}, y);
  }
  else /* prec_order == 0 */
  {
    /* x = s */
    qpb_spinor_xeqy(x, s);
  }

  if(iters == prec_CG_max_iter)
    return -1;

  return iters;
}


int
qpb_preconditioner_CGNE(qpb_spinor_field x, qpb_spinor_field b)
{
  /* Inner CGNE solver invoked by the outer BiCGStab
     (qpb_bicgstab_overlap_kl_pfrac).

     Solves K·x = bprime  in the least-squares sense by iterating the normal
     equations K†K·x = K†·bprime via standard CG, where the operator K and the
     modified RHS bprime are selected by prec_order:

       nPrec=0:  K = M_kernel_op,    bprime = b
                 (CGNE on M·x = b  →  M†M·x = M†·b)

       nPrec=1:  K = M_mult_up_op,   bprime = 3(X^2+1/3)·γ5·b
                 (CGNE on the multiply-up system M''·x = b'  →  M''†M''·x = M''†·b')

     The stopping criterion is on the residual of the normal-equations system
     itself:  ||r||^2 / ||K†·bprime||^2  ≤  prec_CG_epsilon^2 . */

  qpb_spinor_field r      = ov_temp_vecs[8];
  qpb_spinor_field p      = ov_temp_vecs[9];
  qpb_spinor_field w      = ov_temp_vecs[10];
  qpb_spinor_field y      = ov_temp_vecs[11];
  qpb_spinor_field bprime = ov_temp_vecs[12];
  qpb_spinor_field bpp    = ov_temp_vecs[13];

  int iters = 0;

  /* All scalars are real -- K†K is HPD */
  qpb_double bpp_norm, res_norm, new_res_norm, omega, alpha, beta;

  /* --- bprime preprocessing --- */
  if(prec_order == 1)
  {
    /* bprime = 3(X^2+1/3)·γ5·b */
    qpb_spinor_gamma5(w, b);
    X2_shifted_op(y, w, 1.0/3.0);
    qpb_spinor_ax(bprime, (qpb_complex) {3.0, 0.0}, y);
  }
  else /* prec_order == 0 */
  {
    /* bprime = b */
    qpb_spinor_xeqy(bprime, b);
  }

  /* bpp = K†·bprime  (RHS of the normal-equations system K†K·x = K†·bprime) */
  K_inner_conj_op_ptr(bpp, bprime);
  qpb_spinor_xdotx(&bpp_norm, bpp);

  /* --- x0 = 0 --- */
  qpb_spinor_field_set_zero(x);
  qpb_spinor_xeqy(r, bpp);

  /* gamma_0 = ||r0||^2 */
  qpb_spinor_xdotx(&res_norm, r);

  /* p0 = r0 */
  qpb_spinor_xeqy(p, r);

  for(iters = 1; iters < prec_CG_max_iter; iters++)
  {
    /* Stopping on relative residual of the normal-equations system:
       ||r||^2 / ||K†·bprime||^2 <= prec_CG_epsilon^2 */
    if(res_norm / bpp_norm <= prec_CG_epsilon * prec_CG_epsilon)
      break;

    /* Apply K†K to p: w = K·p,  y = K†·w = K†K·p */
    K_inner_op_ptr(w, p);
    K_inner_conj_op_ptr(y, w);

    /* omega = p†·K†K·p = ||K·p||^2 = ||w||^2  (real, positive) */
    qpb_spinor_xdotx(&omega, w);

    /* alpha = ||r||^2 / (p†·K†K·p)  (real, positive) */
    alpha = res_norm / omega;

    /* x += alpha·p */
    qpb_spinor_axpy(x, (qpb_complex) {alpha, 0.}, p, x);

    /* r -= alpha·K†K·p = alpha·y */
    qpb_spinor_axpy(r, (qpb_complex) {-alpha, 0.}, y, r);
    
    /* new_res_norm = ||r_{k+1}||^2 */
    qpb_spinor_xdotx(&new_res_norm, r);

    /* beta = ||r_{k+1}||^2 / ||r_k||^2 */
    beta = new_res_norm / res_norm;

    /* p_{k+1} = r_{k+1} + beta·p_k */
    qpb_spinor_axpy(p, (qpb_complex) {beta, 0.}, p, r);

    res_norm = new_res_norm;
  }

  if(iters == prec_CG_max_iter)
    return -1;

  return iters;
}


/*============================================================================*/
/*                      SECTION 9 - outer solvers                             */
/*============================================================================*/

int
qpb_congrad_overlap_kl_pfrac(qpb_spinor_field x, qpb_spinor_field b,
                                        qpb_double CG_epsilon, int CG_max_iter)
{
  /* Outer CGNE for the overlap operator.
     Solves D_ov·x = b via the normal equations D_ov†D_ov·x = D_ov†·b.

     Preconditioner approximates (D_ov†D_ov)^{-1} via qpb_preconditioner_CG,
     which itself selects the right M (kernel or squared multiply-up) based
     on prec_order.

     If prec_CG_max_iter == 0, no preconditioning is used (s = z). */

  qpb_spinor_field p      = ov_temp_vecs[14];
  qpb_spinor_field r      = ov_temp_vecs[15];
  qpb_spinor_field z      = ov_temp_vecs[16];
  qpb_spinor_field y      = ov_temp_vecs[17];
  qpb_spinor_field w      = ov_temp_vecs[18];
  qpb_spinor_field bprime = ov_temp_vecs[19];
  qpb_spinor_field s      = ov_temp_vecs[20];

  int n_reeval = 100;
  int n_echo = 100;
  int iters = 0;

  qpb_double res_norm, true_res_norm, b_norm, bprime_norm;
  qpb_complex alpha = {1, 0}, omega = {1, 0};
  qpb_complex beta, gamma, new_gamma;

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

  /* Solve M·s0 = z0 for s0 (or copy if no preconditioning) */
  if(prec_CG_max_iter == 0)
    qpb_spinor_xeqy(s, z);
  else
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

    /* Solve M·s = z for s (or copy if no preconditioning) */
    if(prec_CG_max_iter == 0)
      qpb_spinor_xeqy(s, z);
    else
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


int
qpb_bicgstab_overlap_kl_pfrac(qpb_spinor_field x, qpb_spinor_field b,
                               qpb_double epsilon, int max_iter)
{
  /* Outer preconditioned BiCGStab for the overlap operator.
     Solves  D_ov · x = b  directly (no normal equations at the outer level).

     Preconditioner:  K ≈ D_{ov,nPrec} , applied via qpb_preconditioner_CGNE,
     which itself selects the right K (kernel or multiply-up) based on
     prec_order.

     If prec_CG_max_iter == 0, no preconditioning is used (K = I). */

  qpb_spinor_field r0_hat = ov_temp_vecs[14];
  qpb_spinor_field r      = ov_temp_vecs[15];
  qpb_spinor_field p      = ov_temp_vecs[16];
  qpb_spinor_field u      = ov_temp_vecs[17];   /* A·y  or  A·p (unprec) */
  qpb_spinor_field v      = ov_temp_vecs[18];   /* A·z  or  A·s (unprec) */
  qpb_spinor_field y      = ov_temp_vecs[19];   /* K^{-1}·p              */
  qpb_spinor_field z_vec  = ov_temp_vecs[20];   /* K^{-1}·s              */

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
      #if 0  /* LEGACY_BICGSTAB_PREC */
            qpb_preconditioner_CGNE(y, p);
      #else
            qpb_preconditioner_bicgstab(y, p);
      #endif

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
      #if 0  /* LEGACY_BICGSTAB_PREC */
            qpb_preconditioner_CGNE(z_vec, r);
      #else
            qpb_preconditioner_bicgstab(z_vec, r);
      #endif

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


/*============================================================================*/
/*    SECTION 10 - LEGACY / EXPERIMENTAL: inner BiCGStab for nPrec=0          */
/*============================================================================*/
/*
 *  Status: DEAD CODE. Compiled out by `#if 1` below.
 *
 *  Purpose: a one-shot reproduction of the
 *  inner-BiCGStab/outer-BiCGStab combination used in <paper citation>.
 *  Kept in tree so future readers can resurrect the exact numerical
 *  recipe without git archaeology.
 *
 *  Enabling for an experiment requires FOUR edits, both inside `#if 0`
 *  blocks tagged with the marker `LEGACY_BICGSTAB_PREC`:
 *    1. This section: flip `#if 0` → `#if 1` to compile the function.
 *    2. qpb_bicgstab_overlap_kl_pfrac (Section 9): flip the matching
 *       `#if 0` → `#if 1` and the `#if 1` → `#if 0` to swap the inner
 *       call from qpb_preconditioner_CGNE to
 *       qpb_preconditioner_bicgstab.
 *
 *  Constraints when enabled:
 *    - Constraints when enabled: prec_order ∈ {0, 1}. For prec_order=1,
 *      the algorithm uses the multiply-up reorganization (Appendix A.3
 *      of the report) — operator K = M_mult_up_op, bprime =
 *      3(X²+1/3)γ5·b, no post-processing required.
 *
 *  Generated data: <date you ran the experiment, plus a note pointing
 *  to the results file / report section>.
 */

#if 1  /* LEGACY_BICGSTAB_PREC */

static int
qpb_preconditioner_bicgstab(qpb_spinor_field x, qpb_spinor_field b)
{
  if(prec_order != 0 && prec_order != 1)
  {
    error("qpb_preconditioner_bicgstab: prec_order must be 0 or 1, got %d\n",
          prec_order);
    exit(QPB_PARAMETERS_ERROR);
  }

  qpb_spinor_field r      = ov_temp_vecs[8];
  qpb_spinor_field r_hat  = ov_temp_vecs[9];
  qpb_spinor_field p      = ov_temp_vecs[10];
  qpb_spinor_field v      = ov_temp_vecs[11];
  qpb_spinor_field bprime = ov_temp_vecs[12];   /* was unused 's' */
  qpb_spinor_field t      = ov_temp_vecs[13];

  int n_reeval = 100;
  int n_echo = 100;
  int iters = 0;

  qpb_double res_norm, b_norm;
  qpb_complex_double alpha = {1, 0}, omega = {1, 0};
  qpb_complex_double beta, gamma, rho, zeta;
  
  /* --- bprime preprocessing --- */
  if(prec_order == 1)
  {
    /* bprime = 3(X²+1/3)·γ5·b   (matches the CGNE inner's preprocessing
      and the multiply-up derivation in A.3) */
    qpb_spinor_gamma5(v, b);              /* v = γ5·b   (v re-zeroed below) */
    X2_shifted_op(t, v, 1.0/3.0);         /* t = (X²+1/3)·γ5·b              */
    qpb_spinor_ax(bprime, (qpb_complex) {3.0, 0.0}, t);
  }
  else /* prec_order == 0 */
  {
    qpb_spinor_xeqy(bprime, b);
  }

  qpb_double bprime_norm;
  /* ||bprime||^2 */
  qpb_spinor_xdotx(&bprime_norm, bprime);

  /* x0 = 0; zero p,v,t (which were used as scratch above) */
  qpb_spinor_field_set_zero(x);
  qpb_spinor_field_set_zero(p);
  qpb_spinor_field_set_zero(v);
  qpb_spinor_field_set_zero(t);

  /* r0 = bprime - K·x0 = bprime */
  qpb_spinor_xeqy(r, bprime);
  qpb_spinor_xeqy(r_hat, r);

  qpb_spinor_xdotx(&gamma.re, r);
  gamma.im = 0;
  res_norm = gamma.re;
  rho = gamma;

  for(iters = 1; iters < prec_CG_max_iter; iters++)
  {
    if(res_norm / bprime_norm <= prec_CG_epsilon * prec_CG_epsilon)
      break;

    /* gamma = (r̂0, r) */
    qpb_spinor_xdoty(&gamma, r_hat, r);

    /* beta = (gamma/rho)·(alpha/omega) */
    beta = CMUL(CDEV(gamma, rho), CDEV(alpha, omega));

    /* p = r + beta·(p - omega·v) */
    omega = CNEGATE(omega);
    qpb_spinor_axpy(p, omega, v, p);      /* p -= omega·v                */
    qpb_spinor_axpy(p, beta, p, r);       /* p  = beta·p + r             */

    /* v = M_kernel · p */
    K_inner_op_ptr(v, p);

    /* rho = gamma,  alpha = rho / (r̂, v) */
    qpb_spinor_xdoty(&beta, r_hat, v);
    rho = gamma;
    alpha = CDEV(rho, beta);

    /* r -= alpha·v */
    alpha = CNEGATE(alpha);
    qpb_spinor_axpy(r, alpha, v, r);

    /* t = M_kernel · r */
    K_inner_op_ptr(t, r);

    /* omega = (t, r) / (t, t) */
    qpb_spinor_xdoty(&zeta, t, r);
    qpb_spinor_xdotx(&beta.re, t);
    beta.im = 0;
    omega = CDEV(zeta, beta);

    /* x += alpha·y + omega·z  (alpha is currently negated → flip back) */
    alpha = CNEGATE(alpha);
    qpb_spinor_axpy(x, alpha, p, x);     /* x += alpha·p                */
    qpb_spinor_axpy(x, omega, r, x); /* x += omega·r                */

    /* Residual update */
    if(iters % n_reeval == 0)
    {
      K_inner_op_ptr(t, x);
      qpb_spinor_xmy(r, bprime, t);           /* r = b - M_kernel·x              */
    }
    else
    {
      omega = CNEGATE(omega);
      qpb_spinor_axpy(r, omega, t, r);   /* r -= omega·t                */
      omega = CNEGATE(omega);
    }

    qpb_spinor_xdotx(&res_norm, r);
    if(iters % n_echo == 0)
      print(" \t Preconditioner iter = %8d, res = %e\n", iters, res_norm / bprime_norm);
  }

  if(iters == prec_CG_max_iter)
    return -1;

  return iters;
}

#endif /* LEGACY_BICGSTAB_PREC */
