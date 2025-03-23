#ifndef _QPB_BARE_OP_H
#define _QPB_BARE_OP_H 1

#include <qpb_types.h>
#include <qpb_bare_op_defs.h>

void qpb_bare_op_init(void * gauge, qpb_clover_term clover, qpb_double c_sw, \
    qpb_double mass);
void qpb_bare_op_finalize();

void D_bare_op(qpb_spinor_field, qpb_spinor_field);

#endif /* _QPB_BARE_OP_H */