#ifndef _QPB_TWISTED_MASS_LANCZOS_H
#define _QPB_TWISTED_MASS_LANCZOS_H 1
#include <qpb_types.h>

/*
  A bispinor is a flavour doublet of ordinary spinor fields. It owns no
  storage of its own: .flav[0] and .flav[1] are plain qpb_spinor_fields,
  so every existing halo-exchange and dslash routine applies to each
  flavour component unchanged.
*/
typedef struct{
  qpb_spinor_field flav[2];
} qpb_bispinor_field;

qpb_bispinor_field qpb_bispinor_field_init();
void qpb_bispinor_field_finalize(qpb_bispinor_field);

/*
  y = Q_h x, with Q_h the Hermitean non-degenerate twisted mass doublet
  operator (see qpb_twisted_mass_lanczos.c for the derivation).

  y and x must be distinct bispinors; this routine does not work in place.
  qpb_twisted_mass_lanczos_init() must have been called first.
*/
void qpb_tm_doublet_op(qpb_bispinor_field, qpb_bispinor_field, void *,
		       qpb_clover_term, qpb_double, qpb_double,
		       qpb_double, qpb_double);

/*
  Returns |<x,Q_h y> - <Q_h x,y>| / |<x,Q_h y>| for random x, y, which
  vanishes to round-off iff Q_h is Hermitean. Borrows the Lanczos temporary
  vectors, so it must be called after qpb_twisted_mass_lanczos_init() and
  before the first qpb_twisted_mass_lanczos() call.
*/
qpb_double qpb_tm_doublet_hermiticity(void *, qpb_clover_term, qpb_double,
				      qpb_double, qpb_double, qpb_double);

void qpb_twisted_mass_lanczos_init();
/*
  alpha[], beta[], gauge, clover, kappa, c_sw, niter carry the same meaning
  as in qpb_lanczos(); mubar and epsbar are the two extra twisted mass
  parameters. niter > 0 starts a fresh Lanczos run, niter < 0 continues the
  current one by |niter| steps.
*/
void qpb_twisted_mass_lanczos(qpb_double *, qpb_double *, void *,
			      qpb_clover_term, qpb_double, qpb_double,
			      qpb_double, qpb_double, int);
void qpb_twisted_mass_lanczos_finalize();
#endif /* _QPB_TWISTED_MASS_LANCZOS_H */
