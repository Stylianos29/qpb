#ifndef _QPB_OVERLAP_ZOLOTAREV_H
#define _QPB_OVERLAP_ZOLOTAREV_H 1

#include <qpb_types.h>
#include <qpb_kl_defs.h>

void qpb_overlap_Zolotarev_init(void *, qpb_clover_term, int, qpb_double,
     qpb_double, qpb_double, qpb_double,
     qpb_double /*ms_epsilon*/, qpb_double /*prec_ms_epsilon*/, int /*ms_max_iter*/,
     int /*prec_max_iter*/,
     qpb_double /*Lanczos_epsilon*/, int /*Lanczos_max_iters*/,
     qpb_double /*delta_max*/, qpb_double /*delta_min*/);
int  qpb_bicgstab_overlap_Zolotarev(qpb_spinor_field, qpb_spinor_field, qpb_double, int);

// void qpb_overlap_Zolotarev_init(void *, qpb_clover_term, int, qpb_double, qpb_double, qpb_double, qpb_double, qpb_double, int, qpb_double, int, qpb_double, qpb_double);
void qpb_overlap_Zolotarev_finalize();

void qpb_gamma5_sign_function_of_X_Zolotarev(qpb_spinor_field, qpb_spinor_field);
void qpb_overlap_Zolotarev(qpb_spinor_field, qpb_spinor_field);
void qpb_gamma5_overlap_Zolotarev(qpb_spinor_field, qpb_spinor_field);
int qpb_congrad_overlap_Zolotarev(qpb_spinor_field, qpb_spinor_field, qpb_double, int);

#endif /* _QPB_OVERLAP_ZOLOTAREV_H */
