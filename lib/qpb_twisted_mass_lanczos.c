/*
  Lanczos on the square of the Hermitean, non-degenerate twisted mass
  doublet operator. This is the twisted mass counterpart of qpb_lanczos.c
  and deliberately follows its structure step for step.

  ---------------------------------------------------------------------
  The operator
  ---------------------------------------------------------------------

  For a doublet of two quark flavours of equal charge and opposite Wilson
  parameter, the twisted mass Dirac operator in the twisted basis is

      D_h = D_W(m_0) + i*mubar*gamma_5*tau_3 - epsbar*tau_1,

  with tau_i the Pauli matrices in flavour space and D_W the (optionally
  clover improved) Wilson operator as implemented by qpb_apply_dslash() /
  qpb_apply_clover_dslash(), i.e. in the normalisation

      D_W = (m_0 + 4) - (1/2)*hopping - (c_sw/2)*sigma_{mu nu} F_{mu nu},
      m_0 = 1/(2 kappa) - 4.

  The diagonal coefficient is therefore 1/(2 kappa), which is the same
  normalisation used by tmLQCD once its kappa-normalised operator
  (1 + 2 i kappa mu gamma_5 tau_3 - kappa H) is divided by 2 kappa. The
  twisted mass then enters as a bare +i*mu*gamma_5 with no factor of 2 kappa,
  so mubar and epsbar here are directly the a*mu of the ETMC tables.

  The physical (untwisted-basis) masses of the two members of the doublet are

      mu_heavy = mubar + epsbar,     mu_light = mubar - epsbar,

  and positivity of det(D_h) requires mubar > epsbar > 0.

  D_h on its own is not Hermitean, but

      Q_h = gamma_5 * tau_1 * D_h

  is. Writing it out,

      Q_h = (gamma_5 D_W) tau_1 + mubar tau_2 - epsbar gamma_5,

  which is manifestly Hermitean because gamma_5 D_W is (gamma_5-hermiticity
  of the Wilson-clover operator) and tau_1, tau_2, gamma_5 are.

  ---------------------------------------------------------------------
  How it is applied here
  ---------------------------------------------------------------------

  Acting on a bispinor x = (x_0, x_1),

      (Q_h x)_0 = gamma_5 [ (D_W - i*mubar*gamma_5) x_1 - epsbar x_0 ]
      (Q_h x)_1 = gamma_5 [ (D_W + i*mubar*gamma_5) x_0 - epsbar x_1 ].

  Since gamma_5 * (i*mu*gamma_5) = i*mu, pulling gamma_5 through the twisted
  term removes it entirely and leaves

      (Q_h x)_0 = (gamma_5 D_W) x_1 - i*mubar*x_1 - epsbar*gamma_5 x_0
      (Q_h x)_1 = (gamma_5 D_W) x_0 + i*mubar*x_0 - epsbar*gamma_5 x_1.

  This is what qpb_tm_doublet_op() implements. It needs no new kernel: the
  fused gamma_5-dslash wrappers supply (gamma_5 D_W), the twisted term is a
  plain complex axpy, and only the epsbar term costs one qpb_spinor_gamma5()
  per flavour. In particular no assumption is made anywhere about the
  representation of gamma_5 (QPB uses the Dirac basis, in which gamma_5 is
  off-diagonal rather than diag(1,1,-1,-1)).

  ---------------------------------------------------------------------
  What the spectrum of Q_h^2 looks like
  ---------------------------------------------------------------------

  Split D_W = D + B, where D is the hopping term carrying the gamma
  matrices (anti-Hermitean, anticommutes with gamma_5) and B collects the
  Wilson, mass and clover terms (Hermitean, commutes with gamma_5). Setting
  G = mubar*gamma_5*tau_2 + B*tau_1 one has G^2 = mubar^2 + B^2 exactly, and

      Q_h^2 = D^dagger D + (G - epsbar)^2 + [B, D].

  All three terms are Hermitean and the first two are positive semi-definite.
  Since the spectrum of G is +/- sqrt(mubar^2 + b^2), the smallest eigenvalue
  of (G - epsbar)^2 is exactly (mubar - epsbar)^2 = mu_light^2. Hence

      lambda_min(Q_h^2) >= mu_light^2

  holds exactly whenever [B, D] is positive semi-definite, and in particular
  holds exactly in the free theory, where B and D commute. Any violation of
  the bound is carried entirely by the commutator [B, D], which is a lattice
  artefact that vanishes in the continuum. Measuring how close
  lambda_min(Q_h^2) sits to mu_light^2 is the point of this module.

  In the free theory the eigenvalues are known in closed form,

      lambda(p,+/-) = sum_mu sin^2(p_mu)
                      + ( sqrt(mubar^2 + M(p)^2) -/+ epsbar )^2,
      M(p) = m_0 + sum_mu (1 - cos(p_mu)),

  so a run with a unit gauge field, c_sw = 0 and kappa = 1/8 must return
  lambda_min = mu_light^2 to machine precision. That is the cheapest and
  sharpest test of everything in this file.

  ---------------------------------------------------------------------
  Note on precision
  ---------------------------------------------------------------------

  lambda_max is O(50) in this normalisation while lambda_min is O(mu^2),
  which for the physically interesting parameters is O(1e-7). This module is
  meaningless when built with -DSINGLE_PRECISION.
*/

#include <math.h>
#include <qpb_types.h>
#include <qpb_globals.h>
#include <qpb_spinor_field.h>
#include <qpb_spinor_linalg.h>
#include <qpb_comm_halo_spinor_field.h>
#include <qpb_dslash_wrappers.h>
#include <qpb_twisted_mass_lanczos.h>

#define QPB_TM_LANCZOS_NUMB_TEMP_VECS 4

static int n_tm_lanczos;

static qpb_bispinor_field tm_lanczos_temp_vecs[QPB_TM_LANCZOS_NUMB_TEMP_VECS];
static qpb_spinor_field tm_doublet_op_temp_vec;

/*
  Bispinor housekeeping. A bispinor is just a pair of spinor fields, so all
  of these are two calls to the corresponding spinor routine.
*/

qpb_bispinor_field
qpb_bispinor_field_init()
{
  qpb_bispinor_field bispinor_field;
  for(int f=0; f<2; f++)
    bispinor_field.flav[f] = qpb_spinor_field_init();

  return bispinor_field;
}

void
qpb_bispinor_field_finalize(qpb_bispinor_field x)
{
  for(int f=0; f<2; f++)
    qpb_spinor_field_finalize(x.flav[f]);

  return;
}

static void
bispinor_field_set_zero(qpb_bispinor_field x)
{
  for(int f=0; f<2; f++)
    qpb_spinor_field_set_zero(x.flav[f]);

  return;
}

static void
bispinor_field_set_random(qpb_bispinor_field x)
{
  for(int f=0; f<2; f++)
    qpb_spinor_field_set_random(x.flav[f]);

  return;
}

/*
  Bispinor linear algebra. The two inner products each carry out their own
  global reduction; that is two reductions per bispinor dot product instead
  of one, which is negligible next to the four dslash applications that a
  Lanczos step costs.
*/

static void
bispinor_xdotx(qpb_double *dot_prod, qpb_bispinor_field x)
{
  qpb_double d0, d1;
  qpb_spinor_xdotx(&d0, x.flav[0]);
  qpb_spinor_xdotx(&d1, x.flav[1]);
  *dot_prod = d0 + d1;

  return;
}

static void
bispinor_xdoty(qpb_complex_double *dot_prod, qpb_bispinor_field x,
	       qpb_bispinor_field y)
{
  qpb_complex_double d0, d1;
  qpb_spinor_xdoty(&d0, x.flav[0], y.flav[0]);
  qpb_spinor_xdoty(&d1, x.flav[1], y.flav[1]);
  dot_prod->re = d0.re + d1.re;
  dot_prod->im = d0.im + d1.im;

  return;
}

static void
bispinor_ax(qpb_bispinor_field z, qpb_complex alpha, qpb_bispinor_field x)
{
  for(int f=0; f<2; f++)
    qpb_spinor_ax(z.flav[f], alpha, x.flav[f]);

  return;
}

static void
bispinor_axpy(qpb_bispinor_field z, qpb_complex alpha, qpb_bispinor_field x,
	      qpb_bispinor_field y)
{
  for(int f=0; f<2; f++)
    qpb_spinor_axpy(z.flav[f], alpha, x.flav[f], y.flav[f]);

  return;
}

static void
bispinor_xeqy(qpb_bispinor_field x, qpb_bispinor_field y)
{
  for(int f=0; f<2; f++)
    qpb_spinor_xeqy(x.flav[f], y.flav[f]);

  return;
}

/*
  y = Q_h x. See the header comment for the derivation of the three terms.
*/
void
qpb_tm_doublet_op(qpb_bispinor_field y, qpb_bispinor_field x, void *gauge,
		  qpb_clover_term clover, qpb_double mass, qpb_double c_sw,
		  qpb_double mubar, qpb_double epsbar)
{
  void (* gamma5_dslash_func)(qpb_spinor_field, qpb_spinor_field, void **) = NULL;

  void *dslash_args[] =
    {
      gauge,
      &mass,
      &clover,
      &c_sw
    };

  switch(which_dslash_op)
    {
    case QPB_DSLASH_BRILLOUIN:
      if(c_sw)
	gamma5_dslash_func = &qpb_gamma5_clover_bri_dslash;
      else
	gamma5_dslash_func = &qpb_gamma5_bri_dslash;
      break;
    case QPB_DSLASH_STANDARD:
      if(c_sw)
	gamma5_dslash_func = &qpb_gamma5_clover_dslash;
      else
	gamma5_dslash_func = &qpb_gamma5_dslash;
      break;
    }

  /* y_0 = (gamma_5 D_W) x_1,  y_1 = (gamma_5 D_W) x_0 (the tau_1 swap) */
  gamma5_dslash_func(y.flav[0], x.flav[1], dslash_args);
  gamma5_dslash_func(y.flav[1], x.flav[0], dslash_args);

  /* y_0 -= i*mubar*x_1,  y_1 += i*mubar*x_0 (this is mubar*tau_2) */
  qpb_spinor_axpy(y.flav[0], (qpb_complex){0., -mubar}, x.flav[1], y.flav[0]);
  qpb_spinor_axpy(y.flav[1], (qpb_complex){0., +mubar}, x.flav[0], y.flav[1]);

  /* y_f -= epsbar*gamma_5*x_f */
  for(int f=0; f<2; f++)
    {
      qpb_spinor_gamma5(tm_doublet_op_temp_vec, x.flav[f]);
      qpb_spinor_axpy(y.flav[f], (qpb_complex){-epsbar, 0.},
		      tm_doublet_op_temp_vec, y.flav[f]);
    }

  return;
}

qpb_double
qpb_tm_doublet_hermiticity(void *gauge, qpb_clover_term clover,
			   qpb_double kappa, qpb_double c_sw,
			   qpb_double mubar, qpb_double epsbar)
{
  qpb_bispinor_field x = tm_lanczos_temp_vecs[0];
  qpb_bispinor_field y = tm_lanczos_temp_vecs[1];
  qpb_bispinor_field w = tm_lanczos_temp_vecs[2];
  qpb_double mass = 1./(2.*kappa) - 4.;
  qpb_complex_double xqy, qxy;

  bispinor_field_set_random(x);
  bispinor_field_set_random(y);

  qpb_tm_doublet_op(w, y, gauge, clover, mass, c_sw, mubar, epsbar);
  bispinor_xdoty(&xqy, x, w);

  qpb_tm_doublet_op(w, x, gauge, clover, mass, c_sw, mubar, epsbar);
  bispinor_xdoty(&qxy, w, y);

  qpb_double num = sqrt((xqy.re - qxy.re)*(xqy.re - qxy.re) +
			(xqy.im - qxy.im)*(xqy.im - qxy.im));
  qpb_double den = sqrt(xqy.re*xqy.re + xqy.im*xqy.im);

  return (den == 0.) ? num : num/den;
}

void
qpb_twisted_mass_lanczos_init()
{
  for(int i=0; i<QPB_TM_LANCZOS_NUMB_TEMP_VECS; i++)
    {
      tm_lanczos_temp_vecs[i] = qpb_bispinor_field_init();
      bispinor_field_set_zero(tm_lanczos_temp_vecs[i]);
    }
  tm_doublet_op_temp_vec = qpb_spinor_field_init();
  qpb_spinor_field_set_zero(tm_doublet_op_temp_vec);
  qpb_comm_halo_spinor_field_init();

  return;
}

void
qpb_twisted_mass_lanczos_finalize()
{
  for(int i=0; i<QPB_TM_LANCZOS_NUMB_TEMP_VECS; i++)
    {
      qpb_bispinor_field_finalize(tm_lanczos_temp_vecs[i]);
    }
  qpb_spinor_field_finalize(tm_doublet_op_temp_vec);
  qpb_comm_halo_spinor_field_finalize();

  return;
}

/*
  Lanczos on Q_h^2. Q_h is Hermitean, so applying it twice per step gives
  the positive definite Q_h^2 = D_h^dagger D_h directly; this is the same
  trick qpb_lanczos() uses with gamma_5 D_W.

  As in qpb_lanczos(), only the tridiagonal coefficients alpha[] and beta[]
  are accumulated; the Krylov basis is never stored, so there is no
  reorthogonalisation. Loss of orthogonality shows up as repeated copies of
  already converged Ritz values, which does not affect the accuracy of the
  extremal ones (they stay good to O(eps * lambda_max)) and cannot push the
  smallest Ritz value below lambda_min, since Ritz values approach the
  extremes from the inside. Only lambda_min and lambda_max should be read
  off the resulting tridiagonal matrix; the interior values are not
  trustworthy.
*/
void
qpb_twisted_mass_lanczos(qpb_double *alpha, qpb_double *beta, void *gauge,
			 qpb_clover_term clover, qpb_double kappa,
			 qpb_double c_sw, qpb_double mubar, qpb_double epsbar,
			 int niter)
{
  qpb_bispinor_field w = tm_lanczos_temp_vecs[0];
  qpb_bispinor_field r = tm_lanczos_temp_vecs[1];
  qpb_bispinor_field u0 = tm_lanczos_temp_vecs[2];
  qpb_bispinor_field u1 = tm_lanczos_temp_vecs[3];
  qpb_double mass = 1./(2.*kappa) - 4.;

  int restart = (niter>0);
  niter = abs(niter);
  qpb_double a[niter+1], b[niter+1];
  qpb_complex z0;
  qpb_complex_double zd;
  qpb_double norm;

  if(restart)
    {
      bispinor_field_set_random(r);
      bispinor_field_set_zero(u0);
      n_tm_lanczos = 0;

      bispinor_xdotx(&norm, r);
      z0 = (qpb_complex){1./sqrt(norm), 0.};
      bispinor_ax(r, z0, r);
      bispinor_xdotx(&norm, r);
      b[0] = sqrt(norm);
    }
  else
    {
      b[0] = beta[n_tm_lanczos-1];
    }

  for(int i=1; i<=niter; i++)
    {
      z0 = (qpb_complex){1./b[i-1], 0.};
      bispinor_ax(u1, z0, r);
      qpb_tm_doublet_op(w, u1, gauge, clover, mass, c_sw, mubar, epsbar);
      qpb_tm_doublet_op(r, w, gauge, clover, mass, c_sw, mubar, epsbar);

      z0 = (qpb_complex){-b[i-1], 0.};
      bispinor_axpy(r, z0, u0, r);

      bispinor_xdoty(&zd, u1, r);
      a[i] = zd.re;
      z0 = (qpb_complex){-a[i], 0.};
      bispinor_axpy(r, z0, u1, r);

      bispinor_xdotx(&norm, r);
      b[i] = sqrt(norm);
      bispinor_xeqy(u0, u1);
    }

  for(int i=0; i<niter; i++)
    {
      alpha[i+n_tm_lanczos] = a[i+1];
      beta[i+n_tm_lanczos] = b[i+1];
    }

  n_tm_lanczos += niter;

  return;
}
