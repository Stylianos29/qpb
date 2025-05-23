#ifndef _QPB_OVERLAP_KL_PFRAC_H
#define _QPB_OVERLAP_KL_PFRAC_H 1

#include <qpb_types.h>
#include <qpb_kl_defs.h>

void qpb_overlap_kl_pfrac_init(void *, qpb_clover_term, enum qpb_kl_classes, int, qpb_double, qpb_double, qpb_double, qpb_double, qpb_double, int);
void qpb_overlap_kl_pfrac_finalize();

void qpb_gamma5_sign_function_of_X_pfrac(qpb_spinor_field, qpb_spinor_field);
void qpb_overlap_kl_pfrac(qpb_spinor_field, qpb_spinor_field);
void qpb_gamma5_overlap_kl_pfrac(qpb_spinor_field, qpb_spinor_field);
int qpb_congrad_overlap_kl_pfrac(qpb_spinor_field, qpb_spinor_field, qpb_double, int);

#endif /* _QPB_OVERLAP_KL_PFRAC_H */