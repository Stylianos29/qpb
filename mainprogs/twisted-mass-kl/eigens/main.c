/*
  Extremal eigenvalues of Q_h^2, the square of the Hermitean non-degenerate
  twisted mass doublet operator, on a single gauge configuration.

  Only lambda_min and lambda_max are reported. The underlying Lanczos does
  not store the Krylov basis and therefore does not reorthogonalise, so the
  interior Ritz values are contaminated by ghost copies of converged
  eigenvalues and are deliberately not printed. The extremal ones are
  unaffected: they stay accurate to O(eps * lambda_max), which for a double
  precision build is O(1e-14) against a lambda_min of O(1e-7).

  The run converges when the relative change of lambda_min between two
  successive diagonalisations of the Lanczos tridiagonal matrix drops below
  "Lanczos epsilon". Note that Ritz values approach lambda_min from above,
  so an under-converged run overestimates it and never reports a
  spuriously small value.
*/

#include <math.h>
#include <qpb.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>

enum qpb_operators which_dslash_op;

/*
  Eigenvalues of the n-by-n symmetric tridiagonal matrix with diagonal a[]
  and off-diagonal b[], returned sorted in ascending order. Requires n >= 2.

  This costs O(n^3) time and O(n^2) memory per call, which is why it is
  driven by the "Diagonalization frequency" parameter rather than being
  called at every Lanczos iteration.
*/
void
tridiag_eigenv(double *eig, double *a, double *b, int n)
{
  gsl_matrix *A = gsl_matrix_calloc(n, n);
  gsl_matrix_set (A, 0, 0, a[0]);
  gsl_matrix_set (A, 0, 0+1, b[0]);
  for(int i=1; i<n-1; i++)
    {
      gsl_matrix_set(A, i, i, a[i]);
      gsl_matrix_set(A, i, i+1, b[i]);
      gsl_matrix_set(A, i, i-1, b[i-1]);
    }
  gsl_matrix_set(A, n-1, n-1, a[n-1]);
  gsl_matrix_set(A, n-1, n-1-1, b[n-1-1]);

  gsl_vector *e = gsl_vector_alloc(n);
  gsl_eigen_symm_workspace *w = gsl_eigen_symm_alloc(n);
  gsl_eigen_symm(A, e, w);
  gsl_eigen_symm_free(w);
  gsl_matrix_free(A);

  gsl_sort_vector(e);

  for(int i=0; i<n; i++)
    eig[i] = gsl_vector_get(e, i);

  gsl_vector_free(e);
  return;
}

enum {
  CONF_ILDG,
  CONF_RAW_32,
  CONF_RAW_64,
} conf_format;

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

  char aux_string[QPB_MAX_STRING];
  qpb_init_parser(argv[2]);

  int g_dim[ND];
  if(sscanf(qpb_parse("Dimensions"), "%d %d %d %d",
	    g_dim, g_dim+1, g_dim+2, g_dim+3)!=ND)
    {
      error("error parsing for %s\n",
	      "Dimensions");
      exit(QPB_PARSER_ERROR);
    }
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
  strcpy(conf_file, "unit");
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

  /* If Brillouin is selected, parse for whether to project diagonal
     links to SU(N). Default is to not project. In arXiv:1012.3615,
     diagonal links were projected */
  int project_diagonal_links = 0;
  if(which_dslash_op == QPB_DSLASH_BRILLOUIN) {
    char *ret = qpb_parse_optional("Project diagonal links");
    if(ret != NULL) {
      sscanf(ret, "%s", aux_string);
      if(strcmp(aux_string, "yes") == 0)
	project_diagonal_links = 1;
      else if(strcmp(aux_string, "no") == 0)
	project_diagonal_links = 0;
      else {
	error("%s: option should be yes or no\n", "Project diagonal links");
	exit(QPB_PARSER_ERROR);
      }
    } else {
      print(" Will not SU(N)-project diagonal links by default\n");
    }
  }

  qpb_double kappa;
  if(sscanf(qpb_parse("kappa"), "%lf", &kappa)!=1)
    {
      error("error parsing for %s\n",
	    "kappa");
      exit(QPB_PARSER_ERROR);
    }
  qpb_double c_sw;
  if(sscanf(qpb_parse("c_sw"), "%lf", &c_sw)!=1)
    {
      error("error parsing for %s\n",
	    "c_sw");
      exit(QPB_PARSER_ERROR);
    }
  /*
     The two members of the doublet are specified by their own bare twisted
     masses rather than by mubar and epsbar. For the physically interesting
     doublets epsbar/mubar is within a fraction of a percent of one, so
     mu_light = mubar - epsbar suffers a catastrophic cancellation: quoting
     mubar and epsbar to anything less than full precision would corrupt
     exactly the quantity this program exists to measure.
  */
  qpb_double mu_light;
  if(sscanf(qpb_parse("mu_light"), "%lf", &mu_light)!=1)
    {
      error("error parsing for %s\n",
	    "mu_light");
      exit(QPB_PARSER_ERROR);
    }
  qpb_double mu_heavy;
  if(sscanf(qpb_parse("mu_heavy"), "%lf", &mu_heavy)!=1)
    {
      error("error parsing for %s\n",
	    "mu_heavy");
      exit(QPB_PARSER_ERROR);
    }
  qpb_double epsilon;
  if(sscanf(qpb_parse("Lanczos epsilon"), "%lf", &epsilon)!=1)
    {
      error("error parsing for %s\n",
	    "Lanczos epsilon");
      exit(QPB_PARSER_ERROR);
    }
  int max_iters;
  if(sscanf(qpb_parse("Lanczos max iters"), "%d", &max_iters)!=1)
    {
      error("error parsing for %s\n",
	    "Lanczos max iters");
      exit(QPB_PARSER_ERROR);
    }
  /*
     How many Lanczos iterations to run between successive diagonalisations
     of the tridiagonal matrix. Each diagonalisation costs O(n^3) time and
     O(n^2) memory in the current Lanczos dimension n, so at the iteration
     counts this operator needs it must not be done every step. Setting this
     to 1 recovers the behaviour of mainprogs/eigens.
  */
  int diag_freq;
  if(sscanf(qpb_parse("Diagonalization frequency"), "%d", &diag_freq)!=1)
    {
      error("error parsing for %s\n",
	    "Diagonalization frequency");
      exit(QPB_PARSER_ERROR);
    }
  qpb_finalize_parser();

  if(!(mu_heavy > mu_light && mu_light > 0.))
    {
      error("require mu_heavy > mu_light > 0, got mu_light = %g, mu_heavy = %g\n",
	    mu_light, mu_heavy);
      exit(QPB_PARAMETERS_ERROR);
    }
  if(max_iters < 2)
    {
      error("%s should be at least 2\n", "Lanczos max iters");
      exit(QPB_PARAMETERS_ERROR);
    }
  if(diag_freq < 1)
    {
      error("%s should be at least 1\n", "Diagonalization frequency");
      exit(QPB_PARAMETERS_ERROR);
    }

  qpb_double mubar = 0.5*(mu_heavy + mu_light);
  qpb_double epsbar = 0.5*(mu_heavy - mu_light);

  /* initialize cartesian grid and index tables */
  qpb_init(g_dim, procs);

  print(" (Lt, Lz, Ly, Lx) = (%2d,%2d,%2d,%2d)\n",
	problem_params.g_dim[0],
	problem_params.g_dim[1],
	problem_params.g_dim[2],
	problem_params.g_dim[3]);
  print(" Processes = (1,%2d,%2d,%2d)\n", procs[0], procs[1], procs[2]);
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
  print(" kappa = %g\n", kappa);
  print(" m_0 = 1/(2 kappa) - 4 = %+.12e\n", 1./(2.*kappa) - 4.);
  print(" Clover param = %g\n", c_sw);
  print(" mu_light = %+.12e\n", mu_light);
  print(" mu_heavy = %+.12e\n", mu_heavy);
  print(" mubar  = (mu_heavy + mu_light)/2 = %+.12e\n", mubar);
  print(" epsbar = (mu_heavy - mu_light)/2 = %+.12e\n", epsbar);
  print(" Predicted bound mu_light^2 = %+.12e\n", mu_light*mu_light);
  print(" Lanczos epsilon = %g\n", epsilon);
  print(" Lanczos max iters = %d\n", max_iters);
  print(" Diagonalization frequency = %d\n", diag_freq);
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
#ifdef SINGLE_PRECISION
  print(" WARNING: built in single precision. lambda_min is O(mu_light^2)\n");
  print(" WARNING: and cannot be resolved against lambda_max = O(50).\n");
#endif

  qpb_rng_init(seed);
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
  if(ape_niter != 0)
    {
      qpb_gauge_field apegauge = qpb_gauge_field_init();
      print(" APE smear gauge field...\n");
      qpb_apesmear_niter(apegauge, gauge, ape_alpha, ape_niter);
      qpb_gauge_field_copy(gauge, apegauge);
      qpb_gauge_field_finalize(apegauge);

      plaquette = qpb_plaquette(gauge);
      print(" Plaquette = %12.8f\n", plaquette);
    }

  /* Clover term */
  qpb_clover_term clover_term = qpb_clover_term_init();
  qpb_clover_term_get(clover_term, gauge);

  qpb_diagonal_links diagonal_links;
  void *solver_arg_links = NULL;
  switch(which_dslash_op)
    {
    case QPB_DSLASH_BRILLOUIN:
      diagonal_links = qpb_diagonal_links_init();
      qpb_diagonal_links_get(diagonal_links, gauge, project_diagonal_links);
      solver_arg_links = &diagonal_links;
      break;
    case QPB_DSLASH_STANDARD:
      solver_arg_links = &gauge;
      break;
    }

  qpb_double t = qpb_stop_watch(0);
  qpb_twisted_mass_lanczos_init();

  /*
     Q_h is Hermitean by construction. Checking it costs two operator
     applications and validates the flavour structure, the sign of the
     twisted term and the gamma_5 convention all at once, so it is always
     run. It must come before the first qpb_twisted_mass_lanczos() call,
     since it borrows the same temporary vectors.
  */
  qpb_double herm = qpb_tm_doublet_hermiticity(solver_arg_links, clover_term,
					       kappa, c_sw, mubar, epsbar);
  print(" Hermiticity violation of Q_h = %e\n", herm);

  qpb_double *a, *b, *eig;
  a = qpb_alloc(sizeof(qpb_double)*max_iters);
  b = qpb_alloc(sizeof(qpb_double)*max_iters);
  eig = qpb_alloc(sizeof(qpb_double)*max_iters);

  qpb_twisted_mass_lanczos(a, b, solver_arg_links, clover_term, kappa, c_sw,
			   mubar, epsbar, 1);

  qpb_double lambda_min = 0, lambda_max = 0, dlambda = 1, lambda_min0 = 1e3;
  int n = 1, converged = 0;
  for(int i=1; i<max_iters; i++)
    {
      qpb_twisted_mass_lanczos(a, b, solver_arg_links, clover_term, kappa, c_sw,
			       mubar, epsbar, -1);
      n = i+1;
      if(i % diag_freq)
	continue;

      tridiag_eigenv(eig, a, b, n);
      lambda_min = eig[0];
      lambda_max = eig[n-1];
      dlambda = fabs(lambda_min - lambda_min0) / fabs(lambda_min + lambda_min0);
      print(" iter = %6d, lambda_min = %+e, lambda_max = %+e,"
	    " CN = %e (change = %e, target = %e)\n",
	    n, lambda_min, lambda_max, lambda_max/lambda_min, dlambda, epsilon);
      if(dlambda < epsilon*0.5)
	{
	  converged = 1;
	  break;
	}
      lambda_min0 = lambda_min;
    }

  /*
     If the iteration limit was hit part way between two diagonalisations,
     diagonalise once more so that what is reported corresponds to the full
     Krylov space that was actually built.
  */
  if(!converged && (n % diag_freq))
    {
      tridiag_eigenv(eig, a, b, n);
      lambda_min = eig[0];
      lambda_max = eig[n-1];
    }
  t = qpb_stop_watch(t);

  if(!converged)
    {
      print(" WARNING: reached %s = %d without converging lambda_min\n",
	    "Lanczos max iters", max_iters);
    }

  print(" Done, %d Lanczos iterations in t = %g sec\n", n, t);
  print(" lambda_min = %+.12e\n", lambda_min);
  print(" lambda_max = %+.12e\n", lambda_max);
  print(" Condition number = %.12e\n", lambda_max/lambda_min);
  print(" lambda_min / mu_light^2 = %+.12e\n",
	lambda_min/(mu_light*mu_light));
  /*
     One machine readable line per configuration, for the ensemble-wide
     analysis. Fields: conf, kappa, c_sw, mu_light, mu_heavy, lambda_min,
     lambda_max, iterations, converged.
  */
  print(" TM_EIGENS %s %.12e %g %.12e %.12e %.12e %.12e %d %d\n",
	conf_file, kappa, c_sw, mu_light, mu_heavy, lambda_min, lambda_max,
	n, converged);

  free(a);
  free(b);
  free(eig);

  qpb_twisted_mass_lanczos_finalize();
  if(which_dslash_op == QPB_DSLASH_BRILLOUIN)
    qpb_diagonal_links_finalize(diagonal_links);

  qpb_gauge_field_finalize(gauge);
  qpb_clover_term_finalize(clover_term);
  qpb_rng_finalize();
  qpb_finalize();
  return 0;
}
