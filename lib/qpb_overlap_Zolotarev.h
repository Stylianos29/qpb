#ifndef _QPB_OVERLAP_ZOLOTAREV_H
#define _QPB_OVERLAP_ZOLOTAREV_H 1

#include <qpb_types.h>
#include <qpb_kl_defs.h>

/* Rayleigh quotient iteration refinement of the Lanczos estimate of
   lambda_min^2. QPB_RQI_MEASURE refines and reports but leaves the Ritz value
   in use downstream, so the improvement can be measured without perturbing
   anything else in the run. */
enum qpb_RQI_mode {
  QPB_RQI_OFF = 0,
  QPB_RQI_MEASURE = 1,
  QPB_RQI_APPLY = 2
};

void qpb_overlap_Zolotarev_init(void *, qpb_clover_term, int, qpb_double, qpb_double, qpb_double, qpb_double, qpb_double, int, qpb_double, int, qpb_double, qpb_double, int);
void qpb_overlap_Zolotarev_finalize();

void qpb_gamma5_sign_function_of_X_Zolotarev(qpb_spinor_field, qpb_spinor_field);
void qpb_overlap_Zolotarev(qpb_spinor_field, qpb_spinor_field);
void qpb_gamma5_overlap_Zolotarev(qpb_spinor_field, qpb_spinor_field);
int qpb_congrad_overlap_Zolotarev(qpb_spinor_field, qpb_spinor_field, qpb_double, int);

#endif /* _QPB_OVERLAP_ZOLOTAREV_H */
