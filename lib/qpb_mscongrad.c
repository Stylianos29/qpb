#include <math.h>
#include <stdlib.h>
#include <qpb_types.h>
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

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_eigen.h>

#define QPB_MSCONGRAD_NUMB_TEMP_VECS 3
#define MAX_NUMB_SHIFTS 256

/*
 * ============================================================================
 *  Recycled (deflated) multi-shift CG
 * ============================================================================
 *
 *  Standard qpb multi-shift CG augmented with a recycle / deflation subspace U
 *  whose columns approximate the lowest eigenvectors of the square operator
 *  A = gamma5 D gamma5 D = D^dag D, following the "Recycled Multi-Shift Solver"
 *  method note.
 *
 *  U is built once per gauge configuration (lazily, on the first call) by a
 *  Lanczos eigensolver on A; it persists across the subsequent qpb_mscongrad()
 *  calls sharing the same init/finalize bracket (all RHS of one inversion job)
 *  and is freed in qpb_mscongrad_finalize().  Pre-computed once per call:
 *
 *        C ≡ A U ,      M ≡ U^dag C = U^dag A U   (k x k Hermitian),
 *
 *  so that for any shift sigma, H^sigma = U^dag (A + sigma) U = M + sigma I
 *  (U is orthonormalized => U^dag U = I).
 *
 *  TWO changes to the solver, each independently switchable:
 *
 *  (1) DEFLATED INITIAL GUESS  [QPB_RECYCLE_K > 0]:
 *        x^sigma_0 = U (M + sigma I)^{-1} U^dag b   for every shift,
 *        r_0 = b - (A + sigma_0) x^{sigma_0}_0 = b - C z_0 - sigma_0 x_0,
 *        p_0 = p^sigma_0 = r_0 .
 *      (For point sources this is nearly a no-op since U^dag b is tiny; it is
 *       kept because it supplies the range(U) component needed by change (2).)
 *
 *  (2) DEFLATED SEARCH DIRECTION  [QPB_RECYCLE_PROJECT = 1]:
 *        after every base-direction update, project p back into (A'U)^perp:
 *           p <- p - U (H^0)^{-1} U^dag (A' p),   A' = A + sigma_0, H^0 = M + sigma_0 I.
 *      Because A is Hermitian, U^dag(A' p) = (A'U)^dag p = C^dag p + sigma_0 U^dag p,
 *      so NO extra dslash is needed: cost is O(N k) per iteration (a few inner
 *      products, one k x k triangular solve, one rank-k axpy).  This is the
 *      change that actually removes the low modes from the *iteration* (and
 *      hence speeds convergence independently of source/low-mode overlap).
 *
 *  Configurations (override at compile time, e.g. -DQPB_RECYCLE_PROJECT=0):
 *      QPB_RECYCLE_K = 0                      -> plain solver (no recycling).
 *      QPB_RECYCLE_K > 0, QPB_RECYCLE_PROJECT = 0 -> Phase 1 (init-guess only).
 *      QPB_RECYCLE_K > 0, QPB_RECYCLE_PROJECT = 1 -> Phase 2 (full deflated CG).
 *
 *  NOTE on shifted systems: the deflation projector uses the base shift
 *  (sigma_0).  Shifted systems ride along via the unchanged ms-CG zeta
 *  recurrences, so their residuals are no longer exactly collinear with the
 *  base; expect a per-shift residual floor (verify it is below your physics
 *  tolerance).  This is intrinsic to deflating a shared-Krylov multishift.
 *
 *  Tunables:
 *      QPB_RECYCLE_K : number of deflation vectors carried in U.
 *      QPB_RECYCLE_M : Lanczos search dimension used to build U (>= K).
 * ============================================================================
 */
#ifndef QPB_RECYCLE_K
#define QPB_RECYCLE_K 12
#endif
#ifndef QPB_RECYCLE_M
#define QPB_RECYCLE_M 64
#endif
#ifndef QPB_RECYCLE_PROJECT
#define QPB_RECYCLE_PROJECT 1
#endif

#define QPB_RECYCLE_KMAX (QPB_RECYCLE_K > 0 ? QPB_RECYCLE_K : 1)

qpb_spinor_field mscongrad_temp_vecs[QPB_MSCONGRAD_NUMB_TEMP_VECS + MAX_NUMB_SHIFTS];

/*
 * Persistent recycle-subspace state (survives across qpb_mscongrad() calls
 * within one init/finalize bracket).
 */
static qpb_spinor_field recycle_U[QPB_RECYCLE_KMAX];   /* orthonormal U      */
static qpb_spinor_field recycle_C[QPB_RECYCLE_KMAX];   /* C = A U            */
static qpb_complex_double recycle_M[QPB_RECYCLE_KMAX*QPB_RECYCLE_KMAX]; /* M = U^dag C */
static int recycle_k = 0;          /* number of usable columns in U          */
static int recycle_built = 0;      /* has U been built for this config?      */
static int recycle_allocated = 0;  /* are recycle_U/recycle_C allocated?     */

/* ------------------------------------------------------------------------- *
 *  Tiny complex-double helpers for the small (k x k) dense linear algebra.
 * ------------------------------------------------------------------------- */
typedef qpb_complex_double cdbl;

static inline cdbl
cd(double re, double im)
{
  return (cdbl){re, im};
}

static inline cdbl
cd_sub(cdbl a, cdbl b)
{
  return (cdbl){a.re - b.re, a.im - b.im};
}

static inline cdbl
cd_mul(cdbl a, cdbl b)
{
  return (cdbl){a.re*b.re - a.im*b.im, a.re*b.im + a.im*b.re};
}

static inline cdbl
cd_conj(cdbl a)
{
  return (cdbl){a.re, -a.im};
}

static inline cdbl
cd_scale(cdbl a, double s)
{
  return (cdbl){a.re*s, a.im*s};
}

/* convert a complex-double scalar to the solver's working precision */
static inline qpb_complex
to_qc(cdbl a)
{
  return (qpb_complex){(qpb_double)a.re, (qpb_double)a.im};
}

/* ------------------------------------------------------------------------- *
 *  Dense Hermitian Cholesky:  H = M + sigma I = L L^H, L lower-triangular.
 *  (M is k x k Hermitian, positive-definite for sigma > 0, stored row-major.)
 * ------------------------------------------------------------------------- */
static void
recycle_chol_factor(cdbl *L, const cdbl *M, double sigma, int k)
{
  for(int j=0; j<k; j++)
    {
      double d = M[j*k+j].re + sigma;
      for(int q=0; q<j; q++)
	d -= L[j*k+q].re*L[j*k+q].re + L[j*k+q].im*L[j*k+q].im;
      double Ljj = sqrt(d);
      L[j*k+j] = cd(Ljj, 0.);
      for(int i=j+1; i<k; i++)
	{
	  cdbl s = M[i*k+j];           /* off-diagonal: sigma only on diagonal */
	  for(int q=0; q<j; q++)
	    s = cd_sub(s, cd_mul(L[i*k+q], cd_conj(L[j*k+q])));
	  L[i*k+j] = cd_scale(s, 1.0/Ljj);
	}
    }
}

/* Solve L L^H z = t given the Cholesky factor L. */
static void
recycle_chol_solve(cdbl *z, const cdbl *L, const cdbl *t, int k)
{
  cdbl y[QPB_RECYCLE_KMAX];
  for(int i=0; i<k; i++)
    {
      cdbl s = t[i];
      for(int q=0; q<i; q++)
	s = cd_sub(s, cd_mul(L[i*k+q], y[q]));
      y[i] = cd_scale(s, 1.0/L[i*k+i].re);
    }
  for(int i=k-1; i>=0; i--)
    {
      cdbl s = y[i];
      for(int q=i+1; q<k; q++)
	s = cd_sub(s, cd_mul(cd_conj(L[q*k+i]), z[q]));
      z[i] = cd_scale(s, 1.0/L[i*k+i].re);
    }
}

/* Solve (M + sigma I) z = t  (factor + solve in one shot). */
static void
recycle_solve_shift(cdbl *z, const cdbl *M, double sigma, const cdbl *t, int k)
{
  cdbl L[QPB_RECYCLE_KMAX*QPB_RECYCLE_KMAX];
  recycle_chol_factor(L, M, sigma, k);
  recycle_chol_solve(z, L, t, k);
}

/* ------------------------------------------------------------------------- *
 *  Deflated search-direction projection (change (2)):
 *      p <- p - U (H^0)^{-1} U^dag (A' p),
 *  with H^0 = M + sigma0 I prefactored as L0, and
 *      U^dag(A' p) = C^dag p + sigma0 (U^dag p)        (A Hermitian).
 *  No dslash; cost O(N k).
 * ------------------------------------------------------------------------- */
static void
recycle_project(qpb_spinor_field pvec, const cdbl *L0, double sigma0, int k)
{
  cdbl g[QPB_RECYCLE_KMAX], z[QPB_RECYCLE_KMAX];
  for(int a=0; a<k; a++)
    {
      cdbl cdotp, udotp;
      qpb_spinor_xdoty(&cdotp, recycle_C[a], pvec);  /* C_a^dag p */
      qpb_spinor_xdoty(&udotp, recycle_U[a], pvec);  /* U_a^dag p */
      g[a] = cd(cdotp.re + sigma0*udotp.re, cdotp.im + sigma0*udotp.im);
    }
  recycle_chol_solve(z, L0, g, k);
  for(int a=0; a<k; a++)
    qpb_spinor_axpy(pvec, to_qc(cd_scale(z[a], -1.0)), recycle_U[a], pvec);
}

/* ------------------------------------------------------------------------- *
 *  Build the recycle subspace U (and C = A U, M = U^dag C) by Lanczos with
 *  full re-orthogonalization on A = (g5 D)(g5 D), diagonalizing the resulting
 *  tridiagonal, lifting the k lowest Ritz vectors, and orthonormalizing them.
 *  Uses the SAME dslash operator as the solve.
 * ------------------------------------------------------------------------- */
static int
recycle_build(void (*dslash_func)(qpb_spinor_field, qpb_spinor_field, void **),
	      void **dslash_args)
{
  int k = QPB_RECYCLE_K;
  int m = QPB_RECYCLE_M;
  if(k <= 0)
    return 0;
  if(m < k)
    m = k;

  print(" Recycle: building deflation subspace U via Lanczos (k=%d, m=%d)...\n",
	k, m);
  qpb_double tb = qpb_stop_watch(0);

  qpb_spinor_field *lv = qpb_alloc(sizeof(qpb_spinor_field)*m);
  for(int i=0; i<m; i++)
    {
      lv[i] = qpb_spinor_field_init();
      qpb_spinor_field_set_zero(lv[i]);
    }
  qpb_spinor_field av  = qpb_spinor_field_init();
  qpb_spinor_field tmp = qpb_spinor_field_init();
  qpb_spinor_field_set_zero(av);
  qpb_spinor_field_set_zero(tmp);

  double *alpha = qpb_alloc(sizeof(double)*m);
  double *beta  = qpb_alloc(sizeof(double)*m);

  /* v_0 = normalized random vector */
  qpb_spinor_field_set_random(lv[0]);
  qpb_double nrm;
  qpb_spinor_xdotx(&nrm, lv[0]);
  qpb_spinor_ax(lv[0], (qpb_complex){1./sqrt(nrm), 0.}, lv[0]);

  for(int i=0; i<m; i++)
    {
      /* av = A v_i = (g5 D)(g5 D) v_i */
      dslash_func(tmp, lv[i], dslash_args);
      dslash_func(av,  tmp,   dslash_args);

      if(i > 0)
	qpb_spinor_axpy(av, (qpb_complex){-beta[i-1], 0.}, lv[i-1], av);

      qpb_complex_double a;
      qpb_spinor_xdoty(&a, lv[i], av);
      alpha[i] = a.re;
      qpb_spinor_axpy(av, (qpb_complex){-alpha[i], 0.}, lv[i], av);

      /* full re-orthogonalization (twice, for numerical stability) */
      for(int pass=0; pass<2; pass++)
	for(int j=0; j<=i; j++)
	  {
	    qpb_complex_double c;
	    qpb_spinor_xdoty(&c, lv[j], av);
	    qpb_spinor_axpy(av, (qpb_complex){-c.re, -c.im}, lv[j], av);
	  }

      qpb_double bb;
      qpb_spinor_xdotx(&bb, av);
      beta[i] = sqrt(bb);

      if(i < m-1)
	{
	  if(beta[i] < 1e-12)
	    {
	      /* invariant subspace found -- restart with a fresh random
		 vector (it gets re-orthogonalized at the next step) */
	      qpb_spinor_field_set_random(lv[i+1]);
	    }
	  else
	    {
	      qpb_spinor_ax(lv[i+1], (qpb_complex){1./beta[i], 0.}, av);
	    }
	}
    }

  /* diagonalize the symmetric tridiagonal T (alpha[0..m-1], beta[0..m-2]) */
  gsl_matrix *T = gsl_matrix_calloc(m, m);
  for(int i=0; i<m; i++)
    gsl_matrix_set(T, i, i, alpha[i]);
  for(int i=0; i<m-1; i++)
    {
      gsl_matrix_set(T, i, i+1, beta[i]);
      gsl_matrix_set(T, i+1, i, beta[i]);
    }
  gsl_vector *eval = gsl_vector_alloc(m);
  gsl_matrix *evec = gsl_matrix_alloc(m, m);
  gsl_eigen_symmv_workspace *ws = gsl_eigen_symmv_alloc(m);
  gsl_eigen_symmv(T, eval, evec, ws);
  gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_VAL_ASC);
  gsl_eigen_symmv_free(ws);

  /* lift the k lowest Ritz vectors:  U_c = sum_i evec[i][c] v_i */
  for(int c=0; c<k; c++)
    {
      qpb_spinor_field_set_zero(recycle_U[c]);
      for(int i=0; i<m; i++)
	{
	  double sic = gsl_matrix_get(evec, i, c);
	  qpb_spinor_axpy(recycle_U[c], (qpb_complex){sic, 0.}, lv[i],
			  recycle_U[c]);
	}
    }

  print("  lowest Ritz values theta of A = D^dag D :\n");
  for(int c=0; c<k; c++)
    print("   theta[%2d] = %+e\n", c, gsl_vector_get(eval, c));

  gsl_matrix_free(T);
  gsl_vector_free(eval);
  gsl_matrix_free(evec);

  /* modified Gram-Schmidt: enforce U^dag U = I (so H^sigma = M + sigma I) */
  for(int a=0; a<k; a++)
    {
      for(int b=0; b<a; b++)
	{
	  qpb_complex_double c;
	  qpb_spinor_xdoty(&c, recycle_U[b], recycle_U[a]);
	  qpb_spinor_axpy(recycle_U[a], (qpb_complex){-c.re, -c.im},
			  recycle_U[b], recycle_U[a]);
	}
      qpb_double n2;
      qpb_spinor_xdotx(&n2, recycle_U[a]);
      qpb_spinor_ax(recycle_U[a], (qpb_complex){1./sqrt(n2), 0.}, recycle_U[a]);
    }

  /* C = A U */
  for(int c=0; c<k; c++)
    {
      dslash_func(tmp,           recycle_U[c], dslash_args);
      dslash_func(recycle_C[c],  tmp,          dslash_args);
    }

  /* M = U^dag C   (k x k Hermitian) */
  for(int a=0; a<k; a++)
    for(int b=0; b<k; b++)
      {
	qpb_complex_double mm;
	qpb_spinor_xdoty(&mm, recycle_U[a], recycle_C[b]);
	recycle_M[a*k+b] = mm;
      }

  recycle_k = k;
  recycle_built = 1;

  for(int i=0; i<m; i++)
    qpb_spinor_field_finalize(lv[i]);
  free(lv);
  qpb_spinor_field_finalize(av);
  qpb_spinor_field_finalize(tmp);
  free(alpha);
  free(beta);

  tb = qpb_stop_watch(tb);
  print(" Recycle: subspace ready (k=%d), build t = %g secs\n", k, tb);
  return k;
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

  /* allocate the persistent recycle subspace (filled lazily on first solve) */
  if(QPB_RECYCLE_K > 0 && !recycle_allocated)
    {
      for(int i=0; i<QPB_RECYCLE_K; i++)
	{
	  recycle_U[i] = qpb_spinor_field_init();
	  qpb_spinor_field_set_zero(recycle_U[i]);
	  recycle_C[i] = qpb_spinor_field_init();
	  qpb_spinor_field_set_zero(recycle_C[i]);
	}
      recycle_allocated = 1;
      recycle_built = 0;
      recycle_k = 0;
    }
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

  if(recycle_allocated)
    {
      for(int i=0; i<QPB_RECYCLE_K; i++)
	{
	  qpb_spinor_field_finalize(recycle_U[i]);
	  qpb_spinor_field_finalize(recycle_C[i]);
	}
      recycle_allocated = 0;
      recycle_built = 0;
      recycle_k = 0;
    }
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
    p_s[s] = mscongrad_temp_vecs[4+s];

  qpb_double res_norm, b_norm;
  qpb_complex_double alpha_s[ns], alpha;
  qpb_complex_double beta_s[ns], beta0, beta1;
  qpb_complex_double zeta_s[ns][3];
  qpb_complex_double gamma, delta, omega;
  qpb_double mass = 1./(2.*kappa) - 4.;
  void (* dslash_func)(qpb_spinor_field, qpb_spinor_field, void **) = NULL;

  /* deflation state for this solve */
  int do_project = 0;
  cdbl L0[QPB_RECYCLE_KMAX*QPB_RECYCLE_KMAX];   /* Cholesky of H^0 = M + sigma0 I */

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

  /* qpb_spinor_gamma5(r, b); */
  /* dslash_func(p[0], r, dslash_args); */
  /* qpb_spinor_xeqy(b, p[0]); */


  /*
   * Initialize
   */
  qpb_spinor_xdotx(&b_norm, b);

  /*
   * Build the recycle subspace once per configuration (lazy), using the very
   * same dslash operator that the solve uses.
   */
  if(QPB_RECYCLE_K > 0 && !recycle_built)
    recycle_build(dslash_func, dslash_args);

  if(recycle_built && recycle_k > 0)
    {
      /*
       * ---- (1) Deflated initial guess ----
       *
       *   t      = U^dag b
       *   x^j_0  = U (M + sigma_j I)^{-1} t          for every shift j
       *   r_0    = b - (A + sigma_0) x^{sigma_0}_0
       *          = b - C z_0 - sigma_0 x_0           (C = A U, no extra dslash)
       *   p_0    = p^sigma_0 = r_0
       */
      int k = recycle_k;
      qpb_complex_double t[QPB_RECYCLE_KMAX];
      qpb_complex_double z0[QPB_RECYCLE_KMAX];

      for(int a=0; a<k; a++)
	qpb_spinor_xdoty(&t[a], recycle_U[a], b);

      /* base system x[0], physical shift sigmas[0] */
      recycle_solve_shift(z0, recycle_M, sigmas[0], t, k);
      qpb_spinor_field_set_zero(x[0]);
      for(int a=0; a<k; a++)
	qpb_spinor_axpy(x[0], to_qc(z0[a]), recycle_U[a], x[0]);   /* x0 = U z0 */

      /* r0 = b - C z0 - sigma0 x0   (w used as scratch for C z0) */
      qpb_spinor_field_set_zero(w);
      for(int a=0; a<k; a++)
	qpb_spinor_axpy(w, to_qc(z0[a]), recycle_C[a], w);         /* w = C z0  */
      qpb_spinor_xmy(r, b, w);                                    /* r = b - C z0 */
      qpb_spinor_axpy(r, (qpb_complex){-sigmas[0], 0.}, x[0], r); /* r -= s0 x0 */

      /* shifted systems x[s+1], physical shift sigmas[s+1] */
      for(int s=0; s<ns; s++)
	{
	  qpb_complex_double zs[QPB_RECYCLE_KMAX];
	  recycle_solve_shift(zs, recycle_M, sigmas[s+1], t, k);
	  qpb_spinor_field_set_zero(x[s+1]);
	  for(int a=0; a<k; a++)
	    qpb_spinor_axpy(x[s+1], to_qc(zs[a]), recycle_U[a], x[s+1]);
	}

      /* seed all search directions with the (deflated) base residual */
      qpb_spinor_xeqy(p, r);
      for(int s=0; s<ns; s++)
	qpb_spinor_xeqy(p_s[s], r);

      /*
       * ---- (2) Deflated search direction setup ----
       * Prefactor H^0 = M + sigma0 I once, and project the initial base
       * search direction p0 into (A'U)^perp.
       */
      if(QPB_RECYCLE_PROJECT)
	{
	  do_project = 1;
	  recycle_chol_factor(L0, recycle_M, sigmas[0], k);
	  recycle_project(p, L0, sigmas[0], k);
	}
    }
  else
    {
      /* ---- original (no deflation) initialization ---- */
      qpb_spinor_xeqy(r, b);
      qpb_spinor_xeqy(p, b);
      for(int s=0; s<ns; s++)
	qpb_spinor_xeqy(p_s[s], b);
    }

  beta0 = (qpb_complex){1., 0.};
  alpha = (qpb_complex){0., 0.};

  for(int s=0; s<ns; s++)
    {
      zeta_s[s][0] = (qpb_complex){1., 0.};
      zeta_s[s][1] = (qpb_complex){1., 0.};
      alpha_s[s] = (qpb_complex){0., 0.};
    }

  qpb_spinor_xdotx(&res_norm, r);

  if(recycle_built && recycle_k > 0)
    print(" Recycle: deflated initial relative residual = %e, projection = %s\n",
	  res_norm / b_norm, do_project ? "ON" : "OFF");

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
      /*
       * res_norm is the residual of the zero-shift solution
       */
      if(res_norm / b_norm <= epsilon)
	break;

      /*
       * D^+ D on p[0] and add min(shift)
       */
      dslash_count+=1;
      dslash_func(w, p, dslash_args);
      dslash_count+=1;
      dslash_func(y, w, dslash_args);

      qpb_spinor_axpy(w, c_sigma0, p, y);

      qpb_spinor_xdoty(&delta, p, w);
      gamma = (qpb_complex){res_norm, 0.};
      beta1 = CNEGATE(CDEV(gamma, delta));

      /*
       * Compute the zetas, betas and update x
       */
      for(int s=0; s<ns; s++)
	{
	  qpb_complex one = (qpb_complex){1.0, 0.0};
	  qpb_complex c1 = CMUL(zeta_s[s][0], CMUL(zeta_s[s][1], beta0));
	  qpb_complex c2 = CMUL(CMUL(CSUB(zeta_s[s][0], zeta_s[s][1]), beta1), alpha);
	  qpb_complex c3 = CMUL(zeta_s[s][0], CMUL(beta0, (CSUB(one, CMUL(c_sigmas[s], beta1)))));

	  zeta_s[s][2] = CDEV(c1, (CADD(c2, c3)));
    if (CNORM2(zeta_s[s][2]) == 0)
      converged[s] = 1;
	  beta_s[s] = CMUL(beta1, CDEV(zeta_s[s][2], zeta_s[s][1]));

	  /***
	       x[0] is the solution, to the modified operator D^+D + shift0,
	       so there are ns+1 elements in x[], x[1] corresponds to c_sigmas[0], etc.
	  ***/
	  if(!converged[s])
	    qpb_spinor_axpy(x[s+1], CNEGATE(beta_s[s]), p_s[s], x[s+1]);
	}
      qpb_spinor_axpy(x[0], CNEGATE(beta1), p, x[0]);

      /*
       * Update the residual
       */
      if(iters % n_reeval == 0)
	{
	  dslash_count+=1;
    dslash_func(w, x[0], dslash_args);
	  dslash_count+=1;
    dslash_func(y, w, dslash_args);
	  qpb_spinor_axpy(y, c_sigma0, x[0], y);
	  qpb_spinor_xmy(r, b, y);
	}
      else
	{
	  qpb_spinor_axpy(r, beta1, w, r);
	}
      qpb_spinor_xdotx(&omega.re, r);
      omega.im = 0.;
      alpha = CDEV(omega, gamma);
      res_norm = omega.re;

      /*
       * compute alpha_s and update p-vector
       */
      for(int s=0; s<ns; s++)
	{
	  alpha_s[s] = CMUL(alpha, CDEV(CMUL(zeta_s[s][2], beta_s[s]), CMUL(zeta_s[s][1], beta1)));
	  qpb_spinor_axpby(p_s[s], zeta_s[s][2], r, alpha_s[s], p_s[s]);
	}

      qpb_spinor_axpy(p, alpha, p, r);

      /*
       * (2) Deflated search direction: project the updated base p back into
       * (A'U)^perp.  No dslash; cost O(N k).
       */
      if(do_project)
	recycle_project(p, L0, sigmas[0], recycle_k);

      if(iters%n_echo == 0)
	print(" \titers = %8d, res = %e\n", iters, res_norm / b_norm);

      beta0 = beta1;
      for(int s=0; s<ns; s++)
	{
	  zeta_s[s][0] = zeta_s[s][1];
	  zeta_s[s][1] = zeta_s[s][2];
	}

      /*
	Check if one of the shifts has converged
       */
      if(iters % n_check_converged == 0)
	// Check from largest to smallest
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
      error(" residual = %e, relative = %e, t = %g secs\n", res_norm, res_norm / b_norm, t);
      error(" !\n");
      return -1;
    }

  print(" After %d iterations msCG converged, t = %g secs\n", iters, t);
  print(" Total number of dslash applications %d\n", dslash_count);
  print(" Shift = %10g, residual = %e, relative = %e\n", sigmas[0], res_norm, res_norm / b_norm);
  for(int s=0; s<ns; s++)
    print(" Shift = %10g, residual = %e, relative = %e\n", sigmas[s+1], res_s[s], res_s[s] / b_norm);

  if(which_dslash_op == QPB_DSLASH_STANDARD)
    {
      qpb_gauge_field_finalize(gauge_bc);
    }

  return iters;
}
