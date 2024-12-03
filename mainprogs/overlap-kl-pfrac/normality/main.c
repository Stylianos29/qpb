#include <qpb.h>

enum {
  CONF_ILDG,
  CONF_RAW_32,
  CONF_RAW_64,
} conf_format;

enum qpb_operators which_dslash_op;

void
usage(char *argv[])
{
  fprintf(stderr, "Usage: %s geom=NZ,NY,NX PARAM_FILE\n", argv[0]);
  return;
}

int
main(int argc, char *argv[])
{
  /* calls MPI_Init() */
  qpb_msg_passing_init(&argc, &argv);

  /* check and parse command line arguments */
  if(argc != 3)
    {
      usage(argv);
      exit(QPB_ARGS_ERROR);
    }
  char *geom = argv[1];
  int procs[ND-1];
  if(strcmp(strtok(geom, "="), "geom")
     != 0)
    {
      usage(argv);
      exit(QPB_ARGS_ERROR);
    }
  for(int i=0; i<ND-1; i++)
    {
      procs[i] = atoi(strtok(NULL, ","));
    }
  
  /* parse parameter (input) file */
  qpb_init_parser(argv[2]);

  int n_vec;
  if(sscanf(qpb_parse("Number of vectors"), "%d", &n_vec)!=1)
    {
      error("error parsing for %s\n",
	    "Number of vectors");
      exit(QPB_PARSER_ERROR);
    }

  int g_dim[ND];
  if(sscanf(qpb_parse("Dimensions"), "%d %d %d %d",
	    g_dim, g_dim+1, g_dim+2, g_dim+3)!=ND)
    {
      error("error parsing for %s\n", 
	      "Dimensions");
      exit(QPB_PARSER_ERROR);
    }

  char aux_string[QPB_MAX_STRING];
  if(sscanf(qpb_parse("Conf"), "%s", aux_string)!=1)
    {
      error("error parsing for %s\n", 
	    "Conf");
      exit(QPB_PARSER_ERROR);
    }

  enum qpb_field_init_opts conf_opt;
  if(strcmp(aux_string, "file") == 0)
    conf_opt = QPB_FILE;
  else if(strcmp(aux_string, "unit") == 0)
    conf_opt = QPB_UNIT;
  else
    {
      error("%s: option should be one of: ", "Conf");
      error("%s, ", "unit"); 
      error("%s\n", "file"); 
      exit(QPB_PARSER_ERROR);
    };

  char conf_file[QPB_MAX_STRING];
  switch(conf_opt)
    {
    case QPB_ZERO:
      break;
    case QPB_UNIT:
      break;
    case QPB_FILE:
      if(sscanf(qpb_parse("Conf file"), "%s",
		conf_file)!=1)
	{
	  error("error parsing for %s\n", 
		"Conf file");
	  exit(QPB_PARSER_ERROR);
	}  
      if(sscanf(qpb_parse("Conf format"), "%s",
		aux_string)!=1)
	{
	  error("error parsing for %s\n", 
		"Conf format");
	  exit(QPB_PARSER_ERROR);
	}  
      if(strcmp(aux_string, "ildg") == 0)
	conf_format = CONF_ILDG;
      else if(strcmp(aux_string, "raw_32") == 0)
	conf_format = CONF_RAW_32;
      else if(strcmp(aux_string, "raw_64") == 0)
	conf_format = CONF_RAW_64;
      else
	{
	  error("%s: option should be one of: ", "Conf format");
	  error("%s, ", "raw"); 
	  error("%s\n", "ildg"); 
	  exit(QPB_PARSER_ERROR);
	}
      break;
    case QPB_RAND:
      break;
    }

  if(sscanf(qpb_parse("KL class"), "%s", aux_string)!=1)
    {
      error("error parsing for %s\n", 
	    "KL class");
      exit(QPB_PARSER_ERROR);
    }

  enum qpb_kl_classes kl_class;
  if(strcmp(aux_string, "11") == 0)
    kl_class = KL_CLASS_11;
  else
    {
      error("%s: option currently supports: ", "KL class");
      error("%s\n", "11"); 
      exit(QPB_PARSER_ERROR);
    }
  
  int kl_iters;
  if(sscanf(qpb_parse("KL iters"), "%d", &kl_iters)!=1)
    {
      error("error parsing for %s\n", 
	    "KL iters");
      exit(QPB_PARSER_ERROR);
    }
  if(kl_iters<1)
    {
      error("only provide positive integer values for KL iters, quiting\n");
      exit(QPB_PARAMETERS_ERROR);
    }

  qpb_double ape_alpha;
  if(sscanf(qpb_parse("APE alpha"), "%lf", &ape_alpha)!=1)
    {
      error("error parsing for %s\n",
	    "alpha");
      exit(QPB_PARSER_ERROR);
    }

  int ape_niter;
  if(sscanf(qpb_parse("APE iterations"), "%d", &ape_niter)!=1)
    {
      error("error parsing for %s\n",
	    "Smear iterations");
      exit(QPB_PARSER_ERROR);
    }

  int shifts[ND];
  if(sscanf(qpb_parse("Gauge shifts"), "%d %d %d %d",
	    shifts, shifts+1, shifts+2, shifts+3)!=ND)
    {
      error("error parsing for %s\n", 
	      "Shifts");
      exit(QPB_PARSER_ERROR);
    }
  if(shifts[0]<0 ||
     shifts[1]<0 ||
     shifts[2]<0 ||
     shifts[3]<0)
    {
      error("only provide positive shifts, quiting\n");
      exit(QPB_PARAMETERS_ERROR);
    }

  if(shifts[0] > g_dim[0]-1 ||
     shifts[1] > g_dim[1]-1 ||
     shifts[2] > g_dim[2]-1 ||
     shifts[3] > g_dim[3]-1)
    {
      error("shift(s) go beyond lattice length(s), quiting\n");
      exit(QPB_PARAMETERS_ERROR);
    }

  unsigned int seed;
  if(sscanf(qpb_parse("Random seed"), "%u", &seed)!=1)
    {
      error("error parsing for %s\n", 
	    "Random seed");
      exit(QPB_PARSER_ERROR);
    }  

  if(sscanf(qpb_parse("Dslash operator"), "%s", aux_string)!=1)
    {
      error("error parsing for %s\n", 
	    "Dslash operator");
      exit(QPB_PARSER_ERROR);
    }

  /* which_dslash_op is a global */
  if(strcmp(aux_string, "Brillouin") == 0)
    which_dslash_op = QPB_DSLASH_BRILLOUIN;
  else if(strcmp(aux_string, "Standard") == 0)
    which_dslash_op = QPB_DSLASH_STANDARD;
  else
    {
      error("%s: option should be one of: ", "Dslash operator");
      error("%s, ", "Brillouin"); 
      error("%s\n", "Standard"); 
      exit(QPB_PARSER_ERROR);
    };

  qpb_double rho;
  if(sscanf(qpb_parse("rho"), "%lf", &rho)!=1)
    {
      error("error parsing for %s\n",
	    "rho");
      exit(QPB_PARSER_ERROR);
    }

  qpb_double mass;
  if(sscanf(qpb_parse("mass"), "%lf", &mass)!=1)
    {
      error("error parsing for %s\n",
	    "mass");
      exit(QPB_PARSER_ERROR);
    }

  qpb_double c_sw;
  if(sscanf(qpb_parse("c_sw"), "%lf", &c_sw)!=1)
    {
      error("error parsing for %s\n",
	    "c_sw");
      exit(QPB_PARSER_ERROR);
    }

  qpb_double epsilon;
  if(sscanf(qpb_parse("Solver epsilon"), "%lf", &epsilon)!=1)
    {
      error("error parsing for %s\n",
	    "Solver epsilon");
      exit(QPB_PARSER_ERROR);
    }

  int max_iters;
  if(sscanf(qpb_parse("Solver max iters"), "%d", &max_iters)!=1)
    {
      error("error parsing for %s\n",
	    "Solver max iters");
      exit(QPB_PARSER_ERROR);
    }

  qpb_double timebc;
  if(sscanf(qpb_parse("BC in time"), "%lf", &timebc)!=1)
    {
      error("error parsing for %s\n",
	    "BC in time");
      exit(QPB_PARSER_ERROR);
    }
  
  if(timebc != 1 && which_dslash_op == QPB_DSLASH_BRILLOUIN)
    {
      error("\n WARNING: Arbitrary boundary conditions in time are only implemented for Standard Operator\n");
      error(" WARNING: overriding to periodic (setting BC in time = 1)\n\n");
      timebc = 1;
    }

  qpb_double scaling_factor;
  if(sscanf(qpb_parse("Scaling factor"), "%lf", &scaling_factor)!=1)
  {
    error("error parsing for %s\n", "Mu");
    exit(QPB_PARSER_ERROR);
  }

  qpb_finalize_parser();

  /* initialize cartesian grid and index tables */
  qpb_init(g_dim, procs);

  print(" (Lt, Lz, Ly, Lx) = (%2d,%2d,%2d,%2d)\n", 
	problem_params.g_dim[0], 
	problem_params.g_dim[1], 
	problem_params.g_dim[2], 
	problem_params.g_dim[3]);
  print(" Processes = (1,%2d,%2d,%2d)\n", procs[0], procs[1], procs[2]);
  int nthreads = 1;
#ifdef OPENMP
#pragma omp parallel
  nthreads = omp_get_num_threads();
#endif
  print(" Threads per process = %2d\n", nthreads);
  switch(conf_opt)
    {
    case QPB_ZERO:
      print(" Gauge field = Zeros\n");
      break;
    case QPB_UNIT:
      print(" Gauge field = Unit\n");
      break;
    case QPB_FILE:
      if(conf_format == CONF_ILDG)
	{
	  print(" Gauge field (ildg) = %s\n", conf_file);
	}
      else if(conf_format == CONF_RAW_32)
	{
	  print(" Gauge field (raw_32) = %s\n", conf_file);
	}
      else if(conf_format == CONF_RAW_64)
	{
	  print(" Gauge field (raw_64) = %s\n", conf_file);
	}
      break;
    case QPB_RAND:
      print(" Gauge field = Random\n");
      break;
    }

  print(" APE alpha = %g\n", ape_alpha);
  print(" APE iterations = %d\n", ape_niter);
  print(" Conf shifts = %d %d %d %d\n", shifts[0], shifts[1], shifts[2], shifts[3]);
  print(" rho = %g\n", rho);
  print(" Mass = %g\n", mass);
  print(" Clover param = %g\n", c_sw);
  print(" BC in time = %g\n", timebc);
  print(" Solver epsilon = %e\n", epsilon);
  print(" Max solver iters = %d\n", max_iters);
  switch(which_dslash_op)
    {
    case QPB_DSLASH_BRILLOUIN:
      print(" Dslash operator is Brillouin\n");
      break;
    case QPB_DSLASH_STANDARD:
      print(" Dslash operator is Standard\n");
      break;
    }
  qpb_rng_init(seed);
  problem_params.timebc = timebc;

  /* allocate gauge field */
  qpb_gauge_field gauge = qpb_gauge_field_init();

  /* read in configuration */
  switch(conf_opt)
    {
    case QPB_ZERO:
      break;
    case QPB_UNIT:
      qpb_gauge_field_set_unit(gauge);
      break;
    case QPB_FILE:
      if(conf_format == CONF_RAW_32)
	{
	  qpb_read_raw32_gauge(gauge, conf_file);
	}
      else if(conf_format == CONF_RAW_64)
	{
	  qpb_read_raw64_gauge(gauge, conf_file);
	}
      else if(conf_format == CONF_ILDG)
	{
	  qpb_read_ildg_gauge(gauge, conf_file);
	}
      break;
    case QPB_RAND:
      break;
    }

  /* Calculate plaquette */
  qpb_double plaquette = qpb_plaquette(gauge);
  print(" Plaquette = %12.8f\n", plaquette);

  /* APE smear the gauge field */
  qpb_gauge_field apegauge = qpb_gauge_field_init();
  if(ape_niter != 0)
    {
      print(" APE smear gauge field...\n");
      qpb_apesmear_niter(apegauge, gauge, ape_alpha, ape_niter);

      plaquette = qpb_plaquette(apegauge);
      print(" Plaquette = %12.8f\n", plaquette);
    }
  else
    {
      qpb_gauge_field_copy(apegauge, gauge);
    }

  /* Shift it */
  qpb_gauge_field_shift(apegauge, shifts);
  qpb_gauge_field_shift(gauge, shifts);

  /* Clover term */
  qpb_clover_term clover_term = qpb_clover_term_init();
  qpb_clover_term_get(clover_term, apegauge);

  /* Allocate random spinor */
  qpb_spinor_field *eta;
  eta = qpb_alloc(sizeof(qpb_spinor_field)*n_vec);
  for(int i=0; i<n_vec; i++)
    {
      eta[i] = qpb_spinor_field_init();
      qpb_spinor_field_set_gaussian(eta[i]);
    }

  qpb_diagonal_links diagonal_links;
  void *solver_arg_links = NULL;

  switch(which_dslash_op) 
    {
    case QPB_DSLASH_BRILLOUIN:
      diagonal_links = qpb_diagonal_links_init();
      int project_diagonal_links = 0;
      qpb_diagonal_links_get(diagonal_links, apegauge, project_diagonal_links);
      solver_arg_links = &diagonal_links;
      break;
    case QPB_DSLASH_STANDARD:
      solver_arg_links = &apegauge;
      break;
    }
  
  qpb_spinor_field temp_vecs[4];
  for(int i=0; i<4; i++)
    temp_vecs[i] = qpb_spinor_field_init();
  
  qpb_double *diffs;
  diffs = qpb_alloc(sizeof(qpb_double)*n_vec);

  qpb_overlap_kl_pfrac_init(solver_arg_links, clover_term, kl_class, kl_iters, \
                          rho, c_sw, mass, scaling_factor, epsilon, max_iters);

  qpb_double t = qpb_stop_watch(0);
  for(int i=0; i<n_vec; i++)
  {
    print("\n");

    qpb_spinor_field g5Dx = temp_vecs[0];
    qpb_spinor_field g5Dg5Dx = temp_vecs[1];
    /* Compute g5D on eta */
    qpb_overlap_kl_pfrac(g5Dx, eta[i]);
    qpb_spinor_gamma5(g5Dx, g5Dx);
    /* Compute g5D on g5D */
    qpb_overlap_kl_pfrac(g5Dg5Dx, g5Dx);
    qpb_spinor_gamma5(g5Dg5Dx, g5Dg5Dx);

    qpb_spinor_field g5x = temp_vecs[2];
    qpb_spinor_field g5Dg5x = temp_vecs[3];
    qpb_spinor_field Dg5Dg5x = temp_vecs[0];
    /* Compute g5Dg5 on eta */
    qpb_spinor_gamma5(g5x, eta[i]);
    qpb_overlap_kl_pfrac(g5Dg5x, g5x);
    qpb_spinor_gamma5(g5Dg5x, g5Dg5x);
    /* Compute D on g5Dg5 */
    qpb_overlap_kl_pfrac(Dg5Dg5x, g5Dg5x);
    
    qpb_spinor_field x = temp_vecs[2];
    qpb_spinor_xmy(x, Dg5Dg5x, g5Dg5Dx);
    qpb_double x_norm, eta_norm;
    qpb_spinor_xdotx(&x_norm, x);
    qpb_spinor_xdotx(&eta_norm, eta[i]);
    print(" Done vector = %d / %d, ||[D^+D, DD^+]|| = %e\n", i+1, n_vec, x_norm/eta_norm);
    diffs[i] = x_norm/eta_norm;
  }
  t = qpb_stop_watch(t);

  print("\n");
  print(" Done, %d vectors in t = %f sec\n", n_vec, t);
  qpb_overlap_kl_pfrac_finalize();

  print(" ||[D^+D, DD^+]|| (normalized, per stochastic source):\n");
  for(int i=0; i<n_vec; i++)
    print(" %4d %e\n", i, diffs[i]);
  free(diffs);

  for(int i=0; i<4; i++)
    qpb_spinor_field_finalize(temp_vecs[i]);
  
  if(which_dslash_op == QPB_DSLASH_BRILLOUIN)
    qpb_diagonal_links_finalize(diagonal_links);
  
  /* clean up */
  for(int i=0; i<n_vec; i++) 
    {
      qpb_spinor_field_finalize(eta[i]);
    }
  free(eta);
  qpb_gauge_field_finalize(gauge);
  qpb_gauge_field_finalize(apegauge);
  qpb_clover_term_finalize(clover_term);
  qpb_rng_finalize();
  qpb_finalize();
  return 0;
}
