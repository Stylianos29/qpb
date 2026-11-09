/*
  forward: stochastic estimate of the forward application of the Dirac
  operator,

             || D eta ||^2 / || eta ||^2 ,

  for n_vec Gaussian random sources eta, where D is either the standard
  Wilson(-clover) or the Brillouin(-clover) Dirac operator,

             D = D_W + m   (diagonal 4 + m for Standard, 15/8 + m for Brillouin)

  with m the bare mass ("mass" in the parameter file). No matrix sign
  function is involved. The time of the calculation is measured as in
  the other main programs.
*/
#include <qpb.h>
#include <qpb_dslash_wrappers.h>
#include <qpb_comm_halo_spinor_field.h>
#include <qpb_timebc_set_gauge_field.h>

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
	  error("%s, ", "raw_32");
	  error("%s, ", "raw_64");
	  error("%s\n", "ildg");
	  exit(QPB_PARSER_ERROR);
	}
      break;
    case QPB_RAND:
      break;
    }

  qpb_double ape_alpha;
  if(sscanf(qpb_parse("APE alpha"), "%lf", &ape_alpha)!=1)
    {
      error("error parsing for %s\n",
	    "APE alpha");
      exit(QPB_PARSER_ERROR);
    }

  int ape_niter;
  if(sscanf(qpb_parse("APE iterations"), "%d", &ape_niter)!=1)
    {
      error("error parsing for %s\n",
	    "APE iterations");
      exit(QPB_PARSER_ERROR);
    }

  int shifts[ND];
  if(sscanf(qpb_parse("Gauge shifts"), "%d %d %d %d",
	    shifts, shifts+1, shifts+2, shifts+3)!=ND)
    {
      error("error parsing for %s\n",
	      "Gauge shifts");
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

  /* If Brillouin is selected, optionally parse for whether to project
     the diagonal links to SU(N). Default is to not project. */
  int project_diagonal_links = 0;
  if(which_dslash_op == QPB_DSLASH_BRILLOUIN)
    {
      char *ret = qpb_parse_optional("Project diagonal links");
      if(ret != NULL)
	{
	  sscanf(ret, "%s", aux_string);
	  if(strcmp(aux_string, "yes") == 0)
	    project_diagonal_links = 1;
	  else if(strcmp(aux_string, "no") == 0)
	    project_diagonal_links = 0;
	  else
	    {
	      error("%s: option should be yes or no\n", "Project diagonal links");
	      exit(QPB_PARSER_ERROR);
	    }
	}
    }

  /* Bare mass of the Dirac operator: D = D_W + mass */
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
  print(" mass = %g\n", mass);
  print(" Clover param = %g\n", c_sw);
  print(" BC in time = %g\n", timebc);
  switch(which_dslash_op)
  {
  case QPB_DSLASH_BRILLOUIN:
    print(" Dslash operator is Brillouin\n");
    if(project_diagonal_links)
      {
	print(" Will SU(N)-project diagonal links\n");
      }
    else
      {
	print(" Will not SU(N)-project diagonal links\n");
      }
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

  /* Allocate random spinors */
  qpb_spinor_field *eta;
  eta = qpb_alloc(sizeof(qpb_spinor_field)*n_vec);
  for(int i=0; i<n_vec; i++)
    {
      eta[i] = qpb_spinor_field_init();
      qpb_spinor_field_set_gaussian(eta[i]);
    }

  qpb_diagonal_links diagonal_links;
  void *op_arg_links = NULL;

  switch(which_dslash_op)
    {
    case QPB_DSLASH_BRILLOUIN:
      diagonal_links = qpb_diagonal_links_init();
      qpb_diagonal_links_get(diagonal_links, apegauge, project_diagonal_links);
      op_arg_links = &diagonal_links;
      break;
    case QPB_DSLASH_STANDARD:
      op_arg_links = &apegauge;
      break;
    }

  qpb_spinor_field D_eta = qpb_spinor_field_init();

  qpb_double *diffs;
  diffs = qpb_alloc(sizeof(qpb_double)*n_vec);

  /* ---------------------------- overhead ---------------------------- */
  qpb_double t = qpb_stop_watch(0);

  qpb_comm_halo_spinor_field_init();

  /* Set the boundary condition in time on a copy of the gauge field.
     Only implemented for the Standard operator (the Brillouin operator
     has been forced to periodic above) */
  qpb_gauge_field gauge_bc;
  void *gauge_bc_ptr = op_arg_links;
  if(which_dslash_op == QPB_DSLASH_STANDARD)
    {
      gauge_bc = qpb_gauge_field_init();
      qpb_timebc_set_gauge_field(gauge_bc, apegauge, problem_params.timebc);
      gauge_bc_ptr = (void *)&gauge_bc;
    }

  void *dslash_args[] =
    {
      gauge_bc_ptr,
      &mass,
      &clover_term,
      &c_sw
    };

  void (* dslash_func)(qpb_spinor_field, qpb_spinor_field, void **) = NULL;
  switch(which_dslash_op)
    {
    case QPB_DSLASH_BRILLOUIN:
      if(c_sw)
	dslash_func = &qpb_clover_bri_dslash;
      else
	dslash_func = &qpb_bri_dslash;
      break;
    case QPB_DSLASH_STANDARD:
      if(c_sw)
	dslash_func = &qpb_clover_dslash;
      else
	dslash_func = &qpb_dslash;
      break;
    }

  qpb_double t_overhead = qpb_stop_watch(t);
  print(" Total overhead time: %f sec\n", t_overhead);

  /* ------------------------ forward application ------------------------ */
  qpb_double t_forward_sum = 0;
  for(int i=0; i<n_vec; i++)
  {
    print("\n");

    qpb_double t_vec = qpb_stop_watch(0);
    dslash_func(D_eta, eta[i], dslash_args);
    t_vec = qpb_stop_watch(t_vec);
    t_forward_sum += t_vec;

    qpb_double D_eta_norm_squared, eta_norm_squared;
    qpb_spinor_xdotx(&D_eta_norm_squared, D_eta);
    qpb_spinor_xdotx(&eta_norm_squared, eta[i]);

    diffs[i] = D_eta_norm_squared/eta_norm_squared;
    print(" Done vector = %d / %d, ||D(eta)||^2/||eta||^2 = %e, t = %f sec\n",
	  i+1, n_vec, diffs[i], t_vec);
  }
  t = qpb_stop_watch(t);

  print(" Done, %d vectors in t = %f sec\n", n_vec, t);
  print(" Average time per D application = %f sec\n", t_forward_sum/n_vec);

  print(" ||D(eta)||^2/||eta||^2 (per stochastic source):\n");
  for(int i=0; i<n_vec; i++)
    {
      print(" %4d %e\n", i, diffs[i]);
    }
  free(diffs);

  /* clean up */
  qpb_comm_halo_spinor_field_finalize();

  if(which_dslash_op == QPB_DSLASH_STANDARD)
    qpb_gauge_field_finalize(gauge_bc);

  qpb_spinor_field_finalize(D_eta);

  if(which_dslash_op == QPB_DSLASH_BRILLOUIN)
    qpb_diagonal_links_finalize(diagonal_links);

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
