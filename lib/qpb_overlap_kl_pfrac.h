#ifndef _QPB_OVERLAP_KL_PFRAC_H
#define _QPB_OVERLAP_KL_PFRAC_H 1

#include <qpb_types.h>
#include <qpb_kl_defs.h>

void qpb_overlap_kl_pfrac_init(void *, qpb_clover_term, enum qpb_kl_classes,
     int /*kl_iters*/, qpb_double /*rho*/, qpb_double /*c_sw*/,
     qpb_double /*mass*/, qpb_double /*scaling_factor*/,
     qpb_double /*ms_epsilon*/, qpb_double /*prec_ms_epsilon*/, int /*ms_max_iters*/,
     qpb_double /*prec_epsilon*/, int /*prec_max_iter*/,
     qpb_double /*Lanczos_epsilon*/, int /*Lanczos_max_iters*/);
void qpb_overlap_kl_pfrac_finalize();

void qpb_gamma5_sign_function_of_X_pfrac(qpb_spinor_field, qpb_spinor_field);
void qpb_overlap_kl_pfrac(qpb_spinor_field, qpb_spinor_field);
void qpb_gamma5_overlap_kl_pfrac(qpb_spinor_field, qpb_spinor_field);
void qpb_conjugate_overlap_kl_pfrac(qpb_spinor_field, qpb_spinor_field);
int qpb_bicgstab_overlap_kl_pfrac(qpb_spinor_field, qpb_spinor_field, qpb_double, int);

#endif /* _QPB_OVERLAP_KL_PFRAC_H */
