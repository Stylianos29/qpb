#ifndef _QPB_APPLY_DIRAC_LAPLACE_H
#define _QPB_APPLY_DIRAC_LAPLACE_H 1
#include <qpb_types.h>
#include <qpb_spinor_linalg.h>

void qpb_apply_dirac_laplace(qpb_spinor_field, qpb_spinor_field, 
			     qpb_gauge_field, qpb_double);
#endif /* _QPB_APPLY_DIRAC_LAPLACE_H */