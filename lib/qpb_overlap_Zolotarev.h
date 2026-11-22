#ifndef _QPB_OVERLAP_ZOLOTAREV_H
#define _QPB_OVERLAP_ZOLOTAREV_H 1

#include <qpb_types.h>
#include <qpb_kl_defs.h>

/* Hard maxima for the deflation inputs.  A request above any of them is a
   parameter error that stops the program; it is never silently lowered. */
#define QPB_ZOL_DEFL_MAX_K        64  /* k_max of the spectral target       */
#define QPB_ZOL_DEFL_MAX_LAMBDA  0.5  /* lambda_target of the target (|X|)  */

void qpb_overlap_Zolotarev_init(void *, qpb_clover_term, int, qpb_double,
     qpb_double, qpb_double, qpb_double,
     qpb_double /*ms_epsilon*/, qpb_double /*prec_ms_epsilon*/, int /*ms_max_iter*/,
     qpb_double /*prec_epsilon*/, int /*prec_max_iter*/,
     qpb_double /*Lanczos_epsilon*/, int /*Lanczos_max_iters*/,
     qpb_double /*delta_max*/, qpb_double /*delta_min*/);
int  qpb_bicgstab_overlap_Zolotarev(qpb_spinor_field, qpb_spinor_field, qpb_double, int);

// void qpb_overlap_Zolotarev_init(void *, qpb_clover_term, int, qpb_double, qpb_double, qpb_double, qpb_double, qpb_double, int, qpb_double, int, qpb_double, qpb_double);
void qpb_overlap_Zolotarev_finalize();

/* Optional: call before init to choose k from the spectrum instead of the
   compile-time QPB_DEFL_K (smallest k with sqrt(theta_{k+1}) >= lambda_target,
   capped at k_max; lambda_target in |X| units).  Needs
   0 <= lambda_target <= QPB_ZOL_DEFL_MAX_LAMBDA, 0 <= k_max <= QPB_ZOL_DEFL_MAX_K. */
void qpb_overlap_Zolotarev_set_deflation_target(qpb_double /*lambda_target*/,
                                                int /*k_max*/);

/* Optional: call before init to switch deflation off entirely -- no deflation
   build; the Zolotarev interval comes from the plain extreme-eigenvalue
   Lanczos -- whatever QPB_DEFL_K is.  Cannot be combined with a deflation
   target.  Lasts for one init/finalize cycle. */
void qpb_overlap_Zolotarev_disable_deflation(void);

void qpb_gamma5_sign_function_of_X_Zolotarev(qpb_spinor_field, qpb_spinor_field);
void qpb_overlap_Zolotarev(qpb_spinor_field, qpb_spinor_field);
void qpb_gamma5_overlap_Zolotarev(qpb_spinor_field, qpb_spinor_field);
int qpb_congrad_overlap_Zolotarev(qpb_spinor_field, qpb_spinor_field, qpb_double, int);

#endif /* _QPB_OVERLAP_ZOLOTAREV_H */
