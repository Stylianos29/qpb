#include <math.h>
#include <qpb_types.h>
#include <qpb_globals.h>
#include <qpb_spinor_field.h>
#include <qpb_spinor_linalg.h>
#include <qpb_gauge_field.h>
#include <qpb_comm_halo_spinor_field.h>
#include <qpb_comm_halo_gauge_field.h>
#include <qpb_timebc_set_gauge_field.h>
#include <qpb_dslash_wrappers.h>
#include <qpb_stop_watch.h>

#define QPB_MSCONGRAD_NUMB_TEMP_VECS 10   /* was 3; p_s now starts at index 10 */
#define MAX_NUMB_SHIFTS              256
#define MAX_PREC_N                   32   /* hard cap; raise if you want larger N */

qpb_spinor_field mscongrad_temp_vecs[QPB_MSCONGRAD_NUMB_TEMP_VECS + MAX_NUMB_SHIFTS];

/* === Preconditioner config (hard-coded; edit and recompile) === */
static int        prec_N          = 5;    /* Chebyshev order; 0 disables preconditioning */
static qpb_double prec_lambda_min = 0.0;  /* assumed min eigenvalue of X^2 */
static qpb_double prec_lambda_max = 1.63;  /* assumed max eigenvalue of X^2 */

/* === Preconditioner state, recomputed per qpb_mscongrad() call === */
static qpb_double prec_sigma_ref;
static qpb_double prec_alpha;     /* σ_ref + λ_min(X²)  */
static qpb_double prec_beta;      /* σ_ref + λ_max(X²)  */
static qpb_double prec_coeff[MAX_PREC_N];

/* Context captured for the helpers below */
static void (*prec_dslash_func)(qpb_spinor_field, qpb_spinor_field, void **);
static void  *prec_dslash_args_buf[4];


/* Apply (X² + shift)·x → y, using `scratch` for the intermediate γ5·D·x. */
static void
X2_shift_apply(qpb_spinor_field y, qpb_spinor_field x,
               qpb_double shift, qpb_spinor_field scratch)
{
  prec_dslash_func(scratch, x,       (void **)prec_dslash_args_buf);
  prec_dslash_func(y,       scratch, (void **)prec_dslash_args_buf);
  qpb_spinor_axpy(y, (qpb_complex){shift, 0.}, x, y);
}

/* Q = (2(X² + σ_ref) − (β+α)) / (β−α);  maps spectrum [α,β] of (X²+σ_ref) to [-1,1]. */
static void
Q_apply(qpb_spinor_field y, qpb_spinor_field x, qpb_spinor_field scratch)
{
  qpb_double inv_span = 1./(prec_beta - prec_alpha);
  qpb_complex a = { 2.*inv_span, 0. };
  qpb_complex b = { -inv_span*(prec_beta + prec_alpha), 0. };
  X2_shift_apply(y, x, prec_sigma_ref, scratch);   /* y = (X² + σ_ref) x */
  qpb_spinor_axpby(y, a, y, b, x);                 /* y = a·y + b·x       */
}

/* Build c_n for the expansion  1/((β+α)/2 + x(β−α)/2) ≈ Σ c_n T_n(x), x ∈ [-1,1].
   (Same construction as qpb_overlap_Chebyshev.c but with 1/(·) instead of 1/√(·).) */
static void
compute_prec_coefficients(void)
{
  qpb_double xk[MAX_PREC_N], rk[MAX_PREC_N];
  for(int k=0; k<prec_N; k++) {
    xk[k] = cos(M_PI*(k + 0.5)/(qpb_double)prec_N);
    rk[k] = 1./(0.5*(prec_beta + prec_alpha)
               + 0.5*xk[k]*(prec_beta - prec_alpha));
  }
  for(int n=0; n<prec_N; n++) {
    qpb_double sum = 0.;
    for(int k=0; k<prec_N; k++) {
      qpb_double Tn;
      if(n == 0)       Tn = 1.;
      else if(n == 1)  Tn = xk[k];
      else {
        qpb_double Tm2 = 1., Tm1 = xk[k];
        for(int j=2; j<=n; j++) { Tn = 2.*xk[k]*Tm1 - Tm2; Tm2 = Tm1; Tm1 = Tn; }
      }
      sum += rk[k]*Tn;
    }
    qpb_double prefactor = (n == 0 ? 1./prec_N : 2./prec_N);
    prec_coeff[n] = prefactor*sum;
  }
}

/* y ≈ M⁻¹ · x  via Clenshaw on Chebyshev expansion of 1/(X²+σ_ref). */
static void
apply_M_inverse(qpb_spinor_field y, qpb_spinor_field x,
                qpb_spinor_field b1, qpb_spinor_field b2,
                qpb_spinor_field sc)
{
  qpb_complex two = {2., 0.};
  qpb_spinor_field_set_zero(b2);
  qpb_spinor_field_set_zero(b1);
  for(int n = prec_N - 1; n > 0; n--) {
    Q_apply(y, b1, sc);                                          /* y = Q b1   */
    qpb_spinor_ax(y, two, y);                                    /* y = 2 Q b1 */
    qpb_spinor_axpy(y, (qpb_complex){prec_coeff[n], 0.}, x, y);  /* + c_n x    */
    qpb_spinor_xmy(y, y, b2);                                    /* − b2       */
    qpb_spinor_xeqy(b2, b1);
    qpb_spinor_xeqy(b1, y);
  }
  Q_apply(y, b1, sc);
  qpb_spinor_axpy(y, (qpb_complex){prec_coeff[0], 0.}, x, y);
  qpb_spinor_xmy(y, y, b2);
}

void
qpb_mscongrad_init(int numb_shifts)
{
  for(int i=0; i<QPB_MSCONGRAD_NUMB_TEMP_VECS+numb_shifts; i++)
    {
      mscongrad_temp_vecs[i] = qpb_spinor_field_init();
      qpb_spinor_field_set_zero(mscongrad_temp_vecs[i]);
    }
  qpb_comm_halo_spinor_field_init();
  return;
}

void
qpb_mscongrad_finalize(int numb_shifts)
{
  for(int i=0; i<QPB_MSCONGRAD_NUMB_TEMP_VECS+numb_shifts; i++)
    {
      qpb_spinor_field_finalize(mscongrad_temp_vecs[i]);
    }
  qpb_comm_halo_spinor_field_finalize();
  return;
}

/*
 * Original multi-shifted CG tracks solution with zero shift.
 * This version redefines the bare operator to be:  M <- M + s_min, 
 * i.e. adds the smallest shift. Thus we have one less shift, and so
 * we need "n_shift" number of p-vectors instead of "(n_shift+1)".
 */
int
qpb_mscongrad(qpb_spinor_field *x, qpb_spinor_field b, void * gauge,
	      qpb_clover_term clover, qpb_double kappa, int numb_shifts, qpb_double *sigmas,
	      qpb_double c_sw, qpb_double epsilon, int max_iter)
{
  int iters = 0;
  const int n_echo = 10, n_reeval = 100, n_check_converged = 10;

  qpb_spinor_field r = mscongrad_temp_vecs[0];
  qpb_spinor_field y = mscongrad_temp_vecs[1];
  qpb_spinor_field w = mscongrad_temp_vecs[2];
  qpb_spinor_field z       = mscongrad_temp_vecs[4];   /* z   = M r           */
  qpb_spinor_field y_M     = mscongrad_temp_vecs[5];   /* y_M = M p (recursive) */
  qpb_spinor_field v       = mscongrad_temp_vecs[6];   /* v   = A · y_M        */
  qpb_spinor_field cheb_b1 = mscongrad_temp_vecs[7];
  qpb_spinor_field cheb_b2 = mscongrad_temp_vecs[8];
  qpb_spinor_field cheb_sc = mscongrad_temp_vecs[9];

  /*
   *  Sort shifts in ascending order
   */
  for(int i=0; i<numb_shifts; i++){
    for(int j=i; j<numb_shifts; j++){
      if(sigmas[j] < sigmas[i]){
	qpb_double swap = sigmas[i];
	sigmas[i] = sigmas[j];
	sigmas[j] = swap;
      }
    }
  }
  
  // for(int s=0; s<numb_shifts; s++)
  //   {
  //     print(" %g\n", sigmas[s]);
  //   }

  int ns = numb_shifts - 1;

  /*
   *  Offset shifts by smallest shift
   */
  qpb_complex c_sigma0 = (qpb_complex){sigmas[0], 0.};
  qpb_complex c_sigmas[ns];
  for(int s=1; s<numb_shifts; s++)
    c_sigmas[s-1] = (qpb_complex){sigmas[s] - sigmas[0], 0.};

  qpb_spinor_field p_s[ns];
  qpb_spinor_field p;

  p = mscongrad_temp_vecs[3];
  for(int s=0; s<ns; s++)
    p_s[s] = mscongrad_temp_vecs[10+s];
  
  qpb_double res_norm, b_norm;
  qpb_complex_double alpha_s[ns], alpha;
  qpb_complex_double beta_s[ns], beta0, beta1;
  qpb_complex_double zeta_s[ns][3];
  qpb_complex_double gamma, delta, omega;
  qpb_double mass = 1./(2.*kappa) - 4.;
  void (* dslash_func)(qpb_spinor_field, qpb_spinor_field, void **) = NULL;

  /* set boundary condition in time
     !!! currently not implemented for diagonal links !!! */
  void * gauge_bc_ptr;
  qpb_gauge_field gauge_bc;
  if(which_dslash_op == QPB_DSLASH_STANDARD)
    {
      gauge_bc = qpb_gauge_field_init();
      qpb_timebc_set_gauge_field(gauge_bc, *(qpb_gauge_field *)gauge, problem_params.timebc);
      gauge_bc_ptr = (void *)&gauge_bc;
    }
  else
    {
      gauge_bc_ptr = gauge;
    }


  void *dslash_args[] = 
    {
      gauge_bc_ptr,
      &mass,
      &clover,
      &c_sw
    };

  switch(which_dslash_op)
    {
    case QPB_DSLASH_BRILLOUIN:
      if(c_sw)
	dslash_func = &qpb_gamma5_clover_bri_dslash;
      else
	dslash_func = &qpb_gamma5_bri_dslash;	
      break;
    case QPB_DSLASH_STANDARD:
      if(c_sw)
	dslash_func = &qpb_gamma5_clover_dslash;
      else
	dslash_func = &qpb_gamma5_dslash;	
      break;
    }
    
  const int prec_on = (prec_N > 0);

  prec_dslash_args_buf[0] = gauge_bc_ptr;
  prec_dslash_args_buf[1] = &mass;
  prec_dslash_args_buf[2] = &clover;
  prec_dslash_args_buf[3] = &c_sw;
  prec_dslash_func        = dslash_func;

  if(prec_on) {
    prec_sigma_ref = sigmas[numb_shifts - 1];                 /* σ_ref = largest shift */
    print(" \tprec_sigma_ref = %.4e\n", prec_sigma_ref);
    prec_alpha     = prec_sigma_ref + prec_lambda_min;
    prec_beta      = prec_sigma_ref + prec_lambda_max;
    compute_prec_coefficients();
  }

  /* qpb_spinor_gamma5(r, b); */
  /* dslash_func(p[0], r, dslash_args); */
  /* qpb_spinor_xeqy(b, p[0]); */


  /*
   * Initialize
   */
  qpb_double eucl_res_norm;   /* <r,r>   – for the break check          */
  qpb_double res_norm_M;      /* <r,M r> – algorithmic norm (was res_norm) */

  qpb_spinor_xdotx(&b_norm, b);
  qpb_spinor_xeqy(r, b);

  /* z_0 = M r_0,  y_M_0 = M p_0 = M b ; if disabled, M = I. */
  if(prec_on) {
    apply_M_inverse(z, r, cheb_b1, cheb_b2, cheb_sc);
    qpb_spinor_xeqy(y_M, z);
  } else {
    qpb_spinor_xeqy(z, r);
    qpb_spinor_xeqy(y_M, r);
  }

  qpb_spinor_xeqy(p, b);
  for(int s=0; s<ns; s++)
    qpb_spinor_xeqy(p_s[s], b);

  beta0 = (qpb_complex){1., 0.};
  alpha = (qpb_complex){0., 0.};
  for(int s=0; s<ns; s++) {
    zeta_s[s][0] = (qpb_complex){1., 0.};
    zeta_s[s][1] = (qpb_complex){1., 0.};
    alpha_s[s]   = (qpb_complex){0., 0.};
  }

  qpb_spinor_xdotx(&eucl_res_norm, r);
  qpb_spinor_xdoty(&omega, r, z);   omega.im = 0.;   /* ω_0 = <r,M r> = <r,z> */
  res_norm_M = omega.re;

  /* 
     Tracks whether we've converged for each 
     shift independently
  */
  int converged[ns];
  for(int s=0; s<ns; s++)
    converged[s] = 0;
  
  int dslash_count = 0;
  
  qpb_double t = qpb_stop_watch(0);
  for(iters=1; iters<max_iter; iters++)
    {
      /* break on Euclidean residual — unchanged semantics */
      if(eucl_res_norm / b_norm <= epsilon)
        break;

      /* w = A p (unchanged) */
      dslash_count += 1; dslash_func(w, p, dslash_args);
      dslash_count += 1; dslash_func(y, w, dslash_args);
      qpb_spinor_axpy(w, c_sigma0, p, y);

      /* NEW: v = A · y_M  (only needed when preconditioning is on) */
      if(prec_on) {
        dslash_count += 1; dslash_func(v, y_M, dslash_args);
        dslash_count += 1; dslash_func(y, v,   dslash_args);
        qpb_spinor_axpy(v, c_sigma0, y_M, y);
      }

      /* δ = <p, A p>_M = <y_M, w>   (when disabled, y_M == p, so equals <p,w>) */
      qpb_spinor_xdoty(&delta, y_M, w);

      gamma = (qpb_complex){res_norm_M, 0.};   /* γ = ω_k in M-inner product */
      beta1 = CNEGATE(CDEV(gamma, delta));

      /* ζ_s / β_s / x_s / x_0 updates — unchanged */
      for(int s=0; s<ns; s++) {
        qpb_complex one = (qpb_complex){1.0, 0.0};
        qpb_complex c1 = CMUL(zeta_s[s][0], CMUL(zeta_s[s][1], beta0));
        qpb_complex c2 = CMUL(CMUL(CSUB(zeta_s[s][0], zeta_s[s][1]), beta1), alpha);
        qpb_complex c3 = CMUL(zeta_s[s][0], CMUL(beta0, (CSUB(one, CMUL(c_sigmas[s], beta1)))));
        zeta_s[s][2] = CDEV(c1, (CADD(c2, c3)));
        if (CNORM2(zeta_s[s][2]) == 0) converged[s] = 1;
        beta_s[s] = CMUL(beta1, CDEV(zeta_s[s][2], zeta_s[s][1]));
        if(!converged[s])
          qpb_spinor_axpy(x[s+1], CNEGATE(beta_s[s]), p_s[s], x[s+1]);
      }
      qpb_spinor_axpy(x[0], CNEGATE(beta1), p, x[0]);

      /* r and z update — recursive or full refresh */
      if(iters % n_reeval == 0) {
        dslash_count += 1; dslash_func(w, x[0], dslash_args);
        dslash_count += 1; dslash_func(y, w,    dslash_args);
        qpb_spinor_axpy(y, c_sigma0, x[0], y);
        qpb_spinor_xmy(r, b, y);
        if(prec_on) {
          apply_M_inverse(z,   r, cheb_b1, cheb_b2, cheb_sc);   /* refresh z   */
          apply_M_inverse(y_M, p, cheb_b1, cheb_b2, cheb_sc);   /* refresh y_M */
        } else {
          qpb_spinor_xeqy(z,   r);
          qpb_spinor_xeqy(y_M, p);
        }
      } else {
        qpb_spinor_axpy(r, beta1, w, r);                         /* r += β₁ w   */
        if(prec_on) qpb_spinor_axpy(z, beta1, v, z);             /* z += β₁ v   */
        else        qpb_spinor_xeqy(z, r);
      }

      /* Norms */
      qpb_spinor_xdotx(&eucl_res_norm, r);
      qpb_spinor_xdoty(&omega, r, z);   omega.im = 0.;           /* ω_{k+1} = <r, M r> */
      alpha       = CDEV(omega, gamma);
      res_norm_M  = omega.re;

      /* α_s and p_s updates — unchanged */
      for(int s=0; s<ns; s++) {
        alpha_s[s] = CMUL(alpha, CDEV(CMUL(zeta_s[s][2], beta_s[s]), CMUL(zeta_s[s][1], beta1)));
        qpb_spinor_axpby(p_s[s], zeta_s[s][2], r, alpha_s[s], p_s[s]);
      }

      /* p update — unchanged form,  AND recursive y_M = α y_M + z */
      qpb_spinor_axpy(p, alpha, p, r);
      if(prec_on) qpb_spinor_axpy(y_M, alpha, y_M, z);
      else        qpb_spinor_xeqy(y_M, p);

      if(iters%n_echo == 0)
        print(" \titers = %8d, res = %e\n", iters, eucl_res_norm / b_norm);

      beta0 = beta1;
      for(int s=0; s<ns; s++) {
        zeta_s[s][0] = zeta_s[s][1];
        zeta_s[s][1] = zeta_s[s][2];
      }

      /* per-shift convergence check at the end of the iteration — unchanged */
      if(iters % n_check_converged == 0)
        for(int s=ns-1; s>=0; s--)
          if(!converged[s])
	    {
	      qpb_double p_norm;
	      qpb_spinor_xdotx(&p_norm, p_s[s]);
	      // if(CNORM2(beta_s[s])*p_norm < epsilon)
		// {
		  qpb_complex shift = (qpb_complex){sigmas[s+1], 0.};
		  qpb_double res_s;
		  dslash_count+=1;
      dslash_func(y, x[s+1], dslash_args);
		  dslash_count+=1;
      dslash_func(w, y, dslash_args);
		  qpb_spinor_axpy(w, shift, x[s+1], w);
		  qpb_spinor_xmy(y, b, w);
		  qpb_spinor_xdotx(&res_s, y);
		  if(res_s / b_norm <= epsilon)
		    converged[s] = 1;
      else
        break;
		// }	    
	    }
    }
  t = qpb_stop_watch(t);
  qpb_double res_s[ns];

  dslash_count+=1;
  dslash_func(y, x[0], dslash_args);
  dslash_count+=1;
  dslash_func(w, y, dslash_args);
  qpb_spinor_axpy(w, c_sigma0, x[0], w);
  qpb_spinor_xmy(r, b, w);
  qpb_spinor_xdotx(&res_norm, r);
  for(int s=0; s<ns; s++)
    {
      qpb_complex shift = (qpb_complex){sigmas[s+1], 0.};
      dslash_count+=1;
      dslash_func(y, x[s+1], dslash_args);
      dslash_count+=1;
      dslash_func(w, y, dslash_args);
      qpb_spinor_axpy(w, shift, x[s+1], w);
      qpb_spinor_xmy(r, b, w);
      qpb_spinor_xdotx(&res_s[s], r);
    }

  if(iters==max_iter)
    {
      error(" !\n");
      error(" msCG *did not* converge, after %d iterations\n", iters);
      error(" residual = %e, relative = %e, t = %g secs\n", eucl_res_norm, eucl_res_norm / b_norm, t);
      error(" !\n");
      return -1;
    }

  print(" After %d iterations msCG converged, t = %g secs\n", iters, t);
  print(" Total number of dslash applications %d\n", dslash_count);
  print(" Shift = %10g, residual = %e, relative = %e\n", sigmas[0], eucl_res_norm, eucl_res_norm / b_norm);
  for(int s=0; s<ns; s++)
    print(" Shift = %10g, residual = %e, relative = %e\n", sigmas[s+1], res_s[s], res_s[s] / b_norm);
  
  if(which_dslash_op == QPB_DSLASH_STANDARD)
    {
      qpb_gauge_field_finalize(gauge_bc);
    }

  return iters;
}
