#ifndef _QPB_OVERLAP_KL_PFRAC_H
#define _QPB_OVERLAP_KL_PFRAC_H 1

#include <qpb_types.h>
#include <qpb_kl_defs.h>

/* Hard maxima for the deflation inputs.  A request above any of them is a
   parameter error that stops the program; it is never silently lowered. */
#define QPB_KL_DEFL_MAX_VECTORS 1024  /* Lanczos vectors kept in memory     */
#define QPB_KL_DEFL_MAX_K         64  /* k_max of the spectral target       */
#define QPB_KL_DEFL_MAX_LAMBDA   0.5  /* lambda_target of the target (|X|)  */
#define QPB_KL_LM_MAX_K           64  /* modes of the L1 low-mode start     */

void qpb_overlap_kl_pfrac_init(void *, qpb_clover_term, enum qpb_kl_classes,
     int /*kl_iters*/, qpb_double /*rho*/, qpb_double /*c_sw*/,
     qpb_double /*mass*/, qpb_double /*scaling_factor*/,
     qpb_double /*ms_epsilon*/, qpb_double /*prec_ms_epsilon*/, int /*ms_max_iters*/,
     qpb_double /*prec_epsilon*/, int /*prec_max_iter*/,
     qpb_double /*Lanczos_epsilon*/, int /*Lanczos_max_vectors*/);
void qpb_overlap_kl_pfrac_finalize();

/* Optional: call before init to choose k from the spectrum instead of the
   compile-time QPB_DEFL_K (smallest k with sqrt(theta_{k+1}) >= lambda_target,
   capped at k_max; lambda_target in |X| units).  Needs
   0 <= lambda_target <= QPB_KL_DEFL_MAX_LAMBDA, 0 <= k_max <= QPB_KL_DEFL_MAX_K. */
void qpb_overlap_kl_pfrac_set_deflation_target(qpb_double /*lambda_target*/,
                                               int /*k_max*/);

/* Optional: call before init to switch deflation off entirely -- no Lanczos
   build, no extra memory, plain sign function -- whatever QPB_DEFL_K is.
   Cannot be combined with a deflation target.  Lasts for one init/finalize
   cycle. */
void qpb_overlap_kl_pfrac_disable_deflation(void);

/* Optional: call before init to start every L1 (preconditioner) BiCGStab
   from the minimal-residual guess over the k lowest modes of D_p^dag D_p,
   D_p = the L1 overlap (order kl_iters-1, solve mass).  The modes come from
   one Lanczos at init (Ritz-residual target Lanczos_epsilon, at most
   Lanczos_max_vectors vectors).  Needs the preconditioner on,
   1 <= k <= QPB_KL_LM_MAX_K, Lanczos_epsilon > 0 and
   k+1 <= Lanczos_max_vectors <= QPB_KL_DEFL_MAX_VECTORS.  Not calling it
   leaves the solver exactly as before.  Lasts for one init/finalize cycle. */
void qpb_overlap_kl_pfrac_set_lowmode_start(int /*k*/,
                                            qpb_double /*Lanczos_epsilon*/,
                                            int /*Lanczos_max_vectors*/);

void qpb_gamma5_sign_function_of_X_pfrac(qpb_spinor_field, qpb_spinor_field);
void qpb_overlap_kl_pfrac(qpb_spinor_field, qpb_spinor_field);
void qpb_gamma5_overlap_kl_pfrac(qpb_spinor_field, qpb_spinor_field);
void qpb_conjugate_overlap_kl_pfrac(qpb_spinor_field, qpb_spinor_field);
int qpb_bicgstab_overlap_kl_pfrac(qpb_spinor_field, qpb_spinor_field, qpb_double, int);

#endif /* _QPB_OVERLAP_KL_PFRAC_H */
