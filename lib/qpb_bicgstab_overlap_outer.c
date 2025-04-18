#include <math.h>
#include <qpb_types.h>
#include <qpb_globals.h>
#include <qpb_spinor_field.h>
#include <qpb_spinor_linalg.h>
#include <qpb_gauge_field.h>
#include <qpb_comm_halo_spinor_field.h>
#include <qpb_comm_halo_gauge_field.h>
#include <qpb_timebc_set_gauge_field.h>
#include <qpb_dslash_wrappers.h>
#include <qpb_bicgstab.h>
#include <qpb_stop_watch.h>
#include <qpb_kl_defs.h>
#include <qpb_overlap_kl.h>

#define QPB_BICGSTAB_NUMB_TEMP_VECS 7

static qpb_spinor_field bicgstab_temp_vecs[QPB_BICGSTAB_NUMB_TEMP_VECS];
static int precond;
static qpb_overlap_params ov_params;

void
qpb_bicgstab_overlap_outer_init(void * gauge, qpb_clover_term clover, 
				qpb_double rho, qpb_double c_sw, qpb_double mass, int precondition)
{
  if(precondition)
    precond = 1;
  else
    precond = 0;
  
  if(precond)
    qpb_bicgstab_init((char *)"\t", 10000);

  for(int i=0; i<QPB_BICGSTAB_NUMB_TEMP_VECS; i++)
    {
      bicgstab_temp_vecs[i] = qpb_spinor_field_init();
      qpb_spinor_field_set_zero(bicgstab_temp_vecs[i]);
    }
  // qpb_comm_halo_spinor_field_init();
  qpb_overlap_kl_init(gauge, clover, rho, c_sw, mass);
  ov_params.clover = clover;
  ov_params.gauge_ptr = gauge;
  ov_params.rho = rho;
  ov_params.c_sw = c_sw;
  ov_params.mass = mass;
  return;
}

void
qpb_bicgstab_overlap_outer_finalize()
{
  for(int i=0; i<QPB_BICGSTAB_NUMB_TEMP_VECS; i++)
    {
      qpb_spinor_field_finalize(bicgstab_temp_vecs[i]);
    }
  //qpb_comm_halo_spinor_field_finalize();
  qpb_overlap_kl_finalize();
  if(precond)
    qpb_bicgstab_finalize();

  return;
}

#define Op(y, x)							\
  qpb_overlap_kl(y, x, kl_class, kl_iters, epsilon, max_iter);

int
qpb_bicgstab_overlap_outer(qpb_spinor_field x, qpb_spinor_field b,
			   enum qpb_kl_classes kl_class, int kl_iters, qpb_double epsilon, int max_iter)
{
  qpb_spinor_field r0 = bicgstab_temp_vecs[0];
  qpb_spinor_field r = bicgstab_temp_vecs[1];
  qpb_spinor_field p = bicgstab_temp_vecs[2];
  qpb_spinor_field u = bicgstab_temp_vecs[3];
  qpb_spinor_field v = bicgstab_temp_vecs[4];
  // For preconditioned version
  qpb_spinor_field q = bicgstab_temp_vecs[5];
  qpb_spinor_field t = bicgstab_temp_vecs[6];

  int iters = 0;
  const int n_echo = 1, n_reeval = 100;
  qpb_double res_norm, b_norm;
  qpb_complex_double alpha = {1, 0}, omega = {1, 0};
  qpb_complex_double beta, gamma, rho, zeta;
  qpb_spinor_field_set_random(x);
  qpb_spinor_xdotx(&b_norm, b);
  Op(r, x);
  qpb_spinor_xmy(r, b, r);
  qpb_spinor_xeqy(r0, r);
  qpb_spinor_xdotx(&gamma.re, r);
  gamma.im = 0;
  res_norm = gamma.re;
  rho = gamma;
  qpb_double sec = qpb_stop_watch(0);
  for(iters=1; iters<max_iter; iters++)
    {
      if(res_norm / b_norm <= epsilon)
	break;
      
      qpb_spinor_xdoty(&gamma, r0, r);
      beta = CMUL(CDEV(gamma, rho), CDEV(alpha, omega));
      omega = CNEGATE(omega);
      qpb_spinor_axpy(p, omega, u, p);
      qpb_spinor_axpy(p, beta, p, r);
      if(precond) {
	qpb_double mass = ov_params.mass;
	qpb_double rho_ov = ov_params.rho;
	qpb_double factor = 1 - mass / (2*rho_ov);
	qpb_double mass_tilde = mass/factor;
	qpb_double kappa = 1./(8.0+2.0*mass_tilde);
	qpb_bicgstab(q, p, ov_params.gauge_ptr, ov_params.clover, kappa, ov_params.c_sw, 
		     sqrt(epsilon), max_iter);
	Op(u, q);
      }else{
	Op(u, p);
      }
      qpb_spinor_xdoty(&beta, r0, u);
      rho = gamma;
      alpha = CDEV(rho, beta);
      alpha = CNEGATE(alpha);
      qpb_spinor_axpy(r, alpha, u, r);
      if(precond) {
	qpb_double mass = ov_params.mass;
	qpb_double rho_ov = ov_params.rho;
	qpb_double factor = 1 - mass / (2*rho_ov);
	qpb_double mass_tilde = mass/factor;
	qpb_double kappa = 1./(8.0+2.0*mass_tilde);
	qpb_bicgstab(t, r, ov_params.gauge_ptr, ov_params.clover, kappa, ov_params.c_sw, 
		     sqrt(epsilon), max_iter);
	Op(v, t);
      }else{
	Op(v, r);
      }
      qpb_spinor_xdoty(&zeta, v, r);
      qpb_spinor_xdotx(&beta.re, v);
      beta.im = 0;
      omega = CDEV(zeta, beta);
      alpha = CNEGATE(alpha);
      qpb_spinor_axpy(x, omega, r, x);
      qpb_spinor_axpy(x, alpha, p, x);
      if(precond) {       
	qpb_spinor_axpy(x, omega, t, x);
	qpb_spinor_axpy(x, alpha, q, x);		
      }else{
	qpb_spinor_axpy(x, omega, r, x);
	qpb_spinor_axpy(x, alpha, p, x);	
      }
      if(iters % n_reeval == 0)
	{
	  Op(r, x);
	  qpb_spinor_xmy(r, b, r);
	}
      else
	{
	  omega = CNEGATE(omega);
	  qpb_spinor_axpy(r, omega, v, r);
	  omega = CNEGATE(omega);
	}
      qpb_spinor_xdotx(&res_norm, r);
      if((iters % n_echo == 0))
	print(" iters = %8d, res = %e\n", iters, res_norm / b_norm);
    }
  sec = qpb_stop_watch(sec);
  Op(r, x);
  qpb_spinor_xmy(r, b, r);
  qpb_spinor_xdotx(&res_norm, r);
  if(iters==max_iter)
    {
      error(" !\n");
      error(" BiCGStab *did not* converge, after %d iterations\n", iters);
      error(" residual = %e, relative = %e, t = %g secs\n", res_norm, res_norm / b_norm, sec);
      error(" !\n");
      return -1;
    }

  print(" After %d iterations BiCGStab converged\n", iters);
  print(" residual = %e, relative = %e, t = %g secs\n", res_norm, res_norm / b_norm, sec);

  return iters;
}
