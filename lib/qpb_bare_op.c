#include <string.h>
#include <qpb_types.h>
#include <qpb_errors.h>
#include <qpb_globals.h>
#include <qpb_alloc.h>
#include <qpb_spinor_field.h>
#include <qpb_spinor_linalg.h>
#include <qpb_gauge_field.h>
#include <qpb_comm_halo_spinor_field.h>
#include <qpb_comm_halo_gauge_field.h>
#include <qpb_timebc_set_gauge_field.h>
#include <qpb_dslash_wrappers.h>
#include <qpb_stop_watch.h>
#include <qpb_bare_op_defs.h>

static qpb_bare_op_params bare_op_params;


void
qpb_bare_op_init(void * gauge, qpb_clover_term clover, qpb_double c_sw, \
                                                              qpb_double mass)
{
  if(bare_op_params.initialized != QPB_BARE_OP_INITIALIZED)
  {
    qpb_comm_halo_spinor_field_init();

    qpb_gauge_field gauge_bc;
    if(which_dslash_op == QPB_DSLASH_STANDARD)
    {
      gauge_bc = qpb_gauge_field_init();
      qpb_timebc_set_gauge_field(gauge_bc, *(qpb_gauge_field *)gauge,\
                    problem_params.timebc);
      bare_op_params.gauge_ptr = qpb_alloc(sizeof(qpb_gauge_field));
      memcpy(bare_op_params.gauge_ptr, &gauge_bc, sizeof(qpb_gauge_field));
    }
    else
    {
      bare_op_params.gauge_ptr = gauge;
    }

    bare_op_params.c_sw = c_sw;
    bare_op_params.m_bare = mass; // Kernel operator bare mass
    bare_op_params.clover = clover;
    
    switch(which_dslash_op)
    {
    case QPB_DSLASH_BRILLOUIN:
      if(c_sw)
      {
        bare_op_params.g5_dslash_op = &qpb_gamma5_clover_bri_dslash;
        bare_op_params.dslash_op = &qpb_clover_bri_dslash;
      }
      else
      {
        bare_op_params.g5_dslash_op = &qpb_gamma5_bri_dslash;	
        bare_op_params.dslash_op = &qpb_bri_dslash;
      }
      break;
    case QPB_DSLASH_STANDARD:
      if(c_sw)
      {
        bare_op_params.g5_dslash_op = &qpb_gamma5_clover_dslash;
        bare_op_params.dslash_op = &qpb_clover_dslash;
      }
      else
      {
        bare_op_params.g5_dslash_op = &qpb_gamma5_dslash;	
        bare_op_params.dslash_op = &qpb_dslash;	
      }
      break;
    }
    bare_op_params.initialized = QPB_BARE_OP_INITIALIZED;

  }
	
  return;
}


void
qpb_bare_op_finalize()
{
  qpb_comm_halo_spinor_field_finalize();

  if(which_dslash_op == QPB_DSLASH_STANDARD)
    qpb_gauge_field_finalize(*(qpb_gauge_field *)bare_op_params.gauge_ptr);
  
  bare_op_params.initialized = 0;
  
  return;
}


void
D_bare_op(qpb_spinor_field y, qpb_spinor_field x)
{
  /* Implements D - rho*I */

  void *dslash_args[4];

  dslash_args[0] = bare_op_params.gauge_ptr;
  dslash_args[1] = &bare_op_params.m_bare;
  dslash_args[2] = &bare_op_params.clover;
  dslash_args[3] = &bare_op_params.c_sw;
  
  bare_op_params.dslash_op(y, x, dslash_args);
  
  return;
}
