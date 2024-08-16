#ifndef _QPB_OVERLAP_KL_H
#define _QPB_OVERLAP_KL_H 1

#include <qpb_types.h>
#include <qpb_kl_defs.h>

void qpb_overlap_kl_init(void *, qpb_clover_term, qpb_double, qpb_double,\
                                                        qpb_double, qpb_double);
void qpb_overlap_kl_finalize();

void qpb_gamma5_sign_function_partial_fractions(qpb_spinor_field,\
    qpb_spinor_field, qpb_complex, qpb_double*,\
    qpb_double*, int, qpb_double, int);

void qpb_overlap_kl(qpb_spinor_field, qpb_spinor_field, enum qpb_kl_classes,\
                                                        int, qpb_double, int);
void qpb_overlap_kl_single_fraction(qpb_spinor_field, qpb_spinor_field,\
                                    enum qpb_kl_classes, int, qpb_double, int);

void qpb_gamma5_overlap_kl(qpb_spinor_field, qpb_spinor_field,
                                    enum qpb_kl_classes, int, qpb_double, int);

int qpb_congrad_overlap_kl(qpb_spinor_field, qpb_spinor_field,\
                                    enum qpb_kl_classes, int, qpb_double, int);

int qpb_congrad_overlap_kl_partial_fractions(qpb_spinor_field, qpb_spinor_field,\
                                    enum qpb_kl_classes, int, qpb_double, int);

int qpb_inverse_overlap_kl_single_fraction(qpb_spinor_field, qpb_spinor_field,\
                                enum qpb_kl_classes, int, qpb_double, int);


#endif /* _QPB_OVERLAP_KL_H */