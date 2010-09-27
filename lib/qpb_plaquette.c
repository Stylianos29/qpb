#include <qpb_types.h>
#include <qpb_globals.h>
#include <qpb_sun_linalg.h>
#include <qpb_comm_halo_gauge_field.h>

qpb_double 
qpb_plaquette(qpb_gauge_field gauge)
{
  qpb_double plaquette = 0.;
  int lvol = problem_params.l_vol;
  int nprocs = problem_params.nprocs;

  qpb_comm_halo_gauge_field(gauge);

  for(int lv=0; lv<lvol; lv++)
    {
      int v = blk_to_ext[lv];
      for(int mu=0; mu<ND; mu++)
	for(int nu=mu+1; nu<ND; nu++)
	  {
	    qpb_complex *u0, *u1;
	    qpb_complex aux0[NC*NC], aux1[NC*NC];
	    
	    /*
	      u0 = (qpb_complex *)((qpb_link *) gauge.index[v] + mu);
	      u1 = (qpb_complex *)((qpb_link *) gauge.index[nneigh[mu][v]] + nu);
	      U_TIMES_U_00(aux0, u0, u1);
	    
	      u0 = (qpb_complex *)((qpb_link *) gauge.index[nneigh[nu][v]] + mu);
	      U_TIMES_U_01(aux1, aux0, u0);
	    
	      u0 = (qpb_complex *)((qpb_link *) gauge.index[v] + nu);
	      U_TIMES_U_01(aux0, aux1, u0);
	    */
	    
	    /* Plaquette using alternative links. Checks diagonal 
	       neighbors indexing and communication */
	      
	    u0 = (qpb_complex *)((qpb_link *) gauge.index[nneigh[ND+mu][v]] + mu);
	    u1 = (qpb_complex *)((qpb_link *) gauge.index[nneigh[mu+ND][nneigh[nu+ND][v]]] 
				 + nu);
	    U_TIMES_U_11(aux0, u0, u1);
	    
	    u0 = (qpb_complex *)((qpb_link *) gauge.index[nneigh[mu+ND][nneigh[nu+ND][v]]] 
				 + mu);
	    U_TIMES_U_00(aux1, aux0, u0);

	    u0 = (qpb_complex *)((qpb_link *) gauge.index[nneigh[ND+nu][v]] + nu);
	    U_TIMES_U_00(aux0, aux1, u0);

	    plaquette += TRACE_U(aux0);
	  }
    }  
  qpb_double plaquette_vec[nprocs];
  MPI_Gather(&plaquette, sizeof(qpb_double), MPI_BYTE,
	     plaquette_vec, sizeof(qpb_double),
	     MPI_BYTE, QPB_MASTER_PROC, MPI_COMM_WORLD);

  if(am_master)
    {
      plaquette = 0;
      for(int i=0; i<nprocs; i++)
	plaquette += plaquette_vec[i];
    }
  
  plaquette /= ((qpb_double) (nprocs * lvol * ND * (ND - 1))/2);
  MPI_Bcast(&plaquette, sizeof(qpb_double), MPI_BYTE, QPB_MASTER_PROC,
	    MPI_COMM_WORLD);

  return (plaquette);
};