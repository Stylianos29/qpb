#include <qpb.h>
#include <math.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_sort_vector.h>

/*
  Estimates the extreme eigenvalues of H^2 = γ5 Dov,m γ5 Dov,m = Dov,m^† Dov,m,
  with Dov,m the massive overlap operator whose sign function is approximated
  by the partial fraction expansion of the KL (11) class, using the Lanczos
  algorithm (no re-orthogonalization).

  H = γ5 Dov,m is exactly Hermitian also for the approximate sign function,
  up to the precision of the inner multi-shift CG solves. For the exact
  overlap operator: spec(H^2) ⊂ [m^2, 4ρ^2].
*/

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

/*
  Eigenvalues, in ascending order, of the n x n real symmetric tridiagonal
  matrix with diagonal a[0..n-1] and off-diagonal b[0..n-2]
*/
void
tridiag_eigenv(qpb_double *eig, qpb_double *a, qpb_double *b, int n)
{
  gsl_matrix *A = gsl_matrix_calloc(n, n);
  for(int i=0; i<n; i++)
    gsl_matrix_set(A, i, i, a[i]);
  for(int i=0; i<n-1; i++)
    {
      gsl_matrix_set(A, i, i+1, b[i]);
      gsl_matrix_set(A, i+1, i, b[i]);
    }

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

  /* Multi-shift CG (MSCG) parameters, for the inner solves that
     evaluate the sign function. NOTE: the MSCG stopping criterion is on
     the *squared* relative residual. Since Lanczos relies on H being
     Hermitian, a tight value is advisable */
  qpb_double ms_epsilon;
  if(sscanf(qpb_parse("Inner solver epsilon"), "%lf", &ms_epsilon)!=1)
    {
      error("error parsing for %s\n",
	    "Inner solver epsilon");
      exit(QPB_PARSER_ERROR);
    }

  int ms_max_iters;
  if(sscanf(qpb_parse("Inner solver max iters"), "%d", &ms_max_iters)!=1)
    {
      error("error parsing for %s\n",
	    "Inner solver max iters");
      exit(QPB_PARSER_ERROR);
    }

  /* Lanczos parameters */
  qpb_double lanczos_epsilon;
  if(sscanf(qpb_parse("Eigensolver epsilon"), "%lf", &lanczos_epsilon)!=1)
    {
      error("error parsing for %s\n",
	    "Eigensolver epsilon");
      exit(QPB_PARSER_ERROR);
    }

  int lanczos_max_iters;
  if(sscanf(qpb_parse("Eigensolver max iters"), "%d", &lanczos_max_iters)!=1)
    {
      error("error parsing for %s\n",
	    "Eigensolver max iters");
      exit(QPB_PARSER_ERROR);
    }
  if(lanczos_max_iters<2)
    {
      error("Eigensolver max iters should be at least 2, quiting\n");
      exit(QPB_PARAMETERS_ERROR);
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
    error("error parsing for %s\n", "Scaling factor");
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
  print(" mass = %g\n", mass);
  print(" Clover param = %g\n", c_sw);
  print(" BC in time = %g\n", timebc);
  print(" Inner solver epsilon = %e\n", ms_epsilon);
  print(" Inner solver max iters = %d\n", ms_max_iters);
  print(" Eigensolver epsilon = %e\n", lanczos_epsilon);
  print(" Eigensolver max iters = %d\n", lanczos_max_iters);
  switch(which_dslash_op)
    {
    case QPB_DSLASH_BRILLOUIN:
      print(" Dslash operator is Brillouin\n");
      break;
    case QPB_DSLASH_STANDARD:
      print(" Dslash operator is Standard\n");
      break;
    }
  print(" KL iters = %d\n", kl_iters);
  print(" Mu = %.3f\n", scaling_factor);

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

  qpb_overlap_kl_pfrac_init(solver_arg_links, clover_term, KL_CLASS_11, kl_iters, \
                          rho, c_sw, mass, scaling_factor, ms_epsilon, ms_max_iters);

  /* Lanczos vectors */
  qpb_spinor_field r = qpb_spinor_field_init();
  qpb_spinor_field w = qpb_spinor_field_init();
  qpb_spinor_field u0 = qpb_spinor_field_init(); /* previous Lanczos vector */
  qpb_spinor_field u1 = qpb_spinor_field_init(); /* current Lanczos vector */

  /* Diagonal (alpha), off-diagonal (beta) elements of the tridiagonal
     matrix, and its eigenvalues (Ritz values) */
  qpb_double *alpha = qpb_alloc(sizeof(qpb_double)*lanczos_max_iters);
  qpb_double *beta = qpb_alloc(sizeof(qpb_double)*lanczos_max_iters);
  qpb_double *eig = qpb_alloc(sizeof(qpb_double)*lanczos_max_iters);

  /* Random, normalized starting vector */
  qpb_spinor_field_set_random(r);
  qpb_spinor_field_set_zero(u0);
  qpb_double r_norm;
  qpb_spinor_xdotx(&r_norm, r);
  qpb_double beta_prev = sqrt(r_norm);

  qpb_double lambda_min_squared = 0, lambda_max_squared = 0;
  qpb_double lambda_min_squared0 = 0, lambda_max_squared0 = 0;
  qpb_double dlambda_min_squared = 1, dlambda_max_squared = 1;
  int converged = 0;
  int n = 0; /* dimension of the Krylov space built so far */

  qpb_double t = qpb_stop_watch(0);
  for(int k=0; k<lanczos_max_iters; k++)
    {
      /* u1 = r / beta_{k-1} */
      qpb_spinor_ax(u1, (qpb_complex){1./beta_prev, 0.}, r);

      /* r = H^2 u1 = γ5 Dov,m γ5 Dov,m u1 */
      qpb_gamma5_overlap_kl_pfrac(w, u1);
      qpb_gamma5_overlap_kl_pfrac(r, w);

      /* r = r - beta_{k-1} u0 */
      qpb_spinor_axpy(r, (qpb_complex){-beta_prev, 0.}, u0, r);

      /* alpha_k = <u1, r> */
      qpb_complex z;
      qpb_spinor_xdoty(&z, u1, r);
      alpha[k] = z.re;

      /* r = r - alpha_k u1 */
      qpb_spinor_axpy(r, (qpb_complex){-alpha[k], 0.}, u1, r);

      /* beta_k = ||r|| */
      qpb_spinor_xdotx(&r_norm, r);
      beta[k] = sqrt(r_norm);

      qpb_spinor_xeqy(u0, u1);
      beta_prev = beta[k];
      n = k+1;

      /* Ritz values of the (k+1)x(k+1) tridiagonal matrix */
      tridiag_eigenv(eig, alpha, beta, n);
      lambda_min_squared = eig[0];
      lambda_max_squared = eig[n-1];

      if(k > 0)
	{
	  dlambda_min_squared = fabs(lambda_min_squared - lambda_min_squared0)
	    / fabs(lambda_min_squared);
	  dlambda_max_squared = fabs(lambda_max_squared - lambda_max_squared0)
	    / fabs(lambda_max_squared);
	  print(" Lanczos iter = %5d, lambda_min_squared = %e (change = %e), "
		"lambda_max_squared = %e (change = %e), target = %e\n",
		n, lambda_min_squared, dlambda_min_squared,
		lambda_max_squared, dlambda_max_squared,
		lanczos_epsilon);
	  if(dlambda_min_squared < lanczos_epsilon &&
	     dlambda_max_squared < lanczos_epsilon)
	    {
	      converged = 1;
	      break;
	    }
	}

      /* Krylov space became invariant; Ritz values are exact */
      if(beta[k] < 1e-14)
	{
	  print(" Lanczos: invariant subspace reached at iter = %d\n", n);
	  converged = 1;
	  break;
	}

      lambda_min_squared0 = lambda_min_squared;
      lambda_max_squared0 = lambda_max_squared;
    }
  t = qpb_stop_watch(t);

  print("\n");
  if(converged)
    {
      print(" Lanczos converged after %d iters, t = %g sec\n", n, t);
    }
  else
    {
      error(" !\n");
      error(" Lanczos *did not* converge after %d iters, t = %g sec\n", n, t);
      error(" relative changes: lambda_min_squared = %e, "
	    "lambda_max_squared = %e, target = %e\n",
	    dlambda_min_squared, dlambda_max_squared, lanczos_epsilon);
      error(" !\n");
    }
  print(" Extreme eigenvalues of H^2 = (g5 Dov,m)^2, "
	"with exact overlap bounds [m^2, 4 rho^2]:\n");
  print(" lambda_min_squared = %+.12e   (m^2     = %.12e)\n",
	lambda_min_squared, mass*mass);
  print(" lambda_max_squared = %+.12e   (4 rho^2 = %.12e)\n",
	lambda_max_squared, 4*rho*rho);
  print(" Condition number = %e\n", lambda_max_squared/lambda_min_squared);

  free(alpha);
  free(beta);
  free(eig);
  qpb_spinor_field_finalize(r);
  qpb_spinor_field_finalize(w);
  qpb_spinor_field_finalize(u0);
  qpb_spinor_field_finalize(u1);

  qpb_overlap_kl_pfrac_finalize();

  if(which_dslash_op == QPB_DSLASH_BRILLOUIN)
    qpb_diagonal_links_finalize(diagonal_links);

  /* clean up */
  qpb_gauge_field_finalize(gauge);
  qpb_gauge_field_finalize(apegauge);
  qpb_clover_term_finalize(clover_term);
  qpb_rng_finalize();
  qpb_finalize();
  return 0;
}
