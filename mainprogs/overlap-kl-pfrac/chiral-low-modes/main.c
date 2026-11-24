#include <qpb.h>
#include <math.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_sort_vector.h>

/*
  Phase A of the deflated-preconditioner plan: the low modes of the overlap
  operator from its chiral structure, for one configuration.

  With s = sign_n(X) the KL (11) partial fraction at order n = "KL iters",
  P± = (1 ± g5)/2 and D = a + b g5 s, a = rho + am/2, b = rho - am/2:

     D^dag D = (1 - (am)^2/(4 rho^2)) K + (am)^2 - b^2 (1 - s^2),
     K|±     = 2 rho^2 (1 + A±),        A± = P± g5 s P±,

  so the lowest eigenvectors of A± (eigenvalues mu near -1) are, up to
  O(1 - s^2), the low modes of D for EVERY mass.  This program:
    1. checks that s is Hermitian (A± must be, for Lanczos to be valid);
    2. runs a Lanczos with full re-orthogonalization on A+ and on A-,
       keeping the "Chiral modes" lowest Ritz pairs of each sector;
    3. for every am in "Masses", measures lambda^2 = ||D u||^2 and the
       residual ||D^dag D u - lambda^2 u|| of each kept mode u.

  The sign function is built as in the solvers ("Lanczos epsilon",
  "Lanczos max vectors" and "Deflation = on/off" drive its own deflation);
  with deflation on, s is Hermitian only up to the Ritz residuals of that
  build, which must therefore converge.  The chiral Lanczos has its own keys.
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

/* Hard limits of the chiral-mode inputs (stop on violation, never clamp). */
#define CHIRAL_MAX_MODES  64    /* "Chiral modes" per sector               */
#define CHIRAL_MAX_MASSES 64    /* entries of "Masses"                     */


/* s x = g5 (g5 s x):  the library applies g5 sign(X). */
static void
apply_sign(qpb_spinor_field y, qpb_spinor_field x)
{
  qpb_gamma5_sign_function_of_X_pfrac(y, x);
  qpb_spinor_gamma5(y, y);
}


/* y = P_sigma x = (x + sigma g5 x) / 2,  sigma = +1 or -1.  y may alias x;
   tmp must differ from both. */
static void
chiral_project(qpb_spinor_field y, qpb_spinor_field x, int sigma,
	       qpb_spinor_field tmp)
{
  qpb_spinor_gamma5(tmp, x);
  qpb_spinor_axpby(y, (qpb_complex){0.5, 0.}, x,
		   (qpb_complex){0.5*sigma, 0.}, tmp);
}


/* y = A_sigma x = P_sigma g5 s P_sigma x  (one sign-function application).
   t1, t2 are scratch, distinct from x and y. */
static void
apply_chiral(qpb_spinor_field y, qpb_spinor_field x, int sigma,
	     qpb_spinor_field t1, qpb_spinor_field t2)
{
  chiral_project(t1, x, sigma, t2);              /* t1 = P x          */
  qpb_gamma5_sign_function_of_X_pfrac(t2, t1);   /* t2 = g5 s P x     */
  chiral_project(y, t2, sigma, t1);              /* y  = P g5 s P x   */
}


/*
  Lanczos with full re-orthogonalization on A_sigma, restricted to the chiral
  sector P_sigma.  Stops when the Ritz-residual bound beta_m |S_{m,i}| of
  each of the nmodes lowest Ritz pairs is below eps (absolute: ||A|| ~ 1), or
  at mmax vectors.  The nmodes lowest Ritz vectors are lifted into u[] (which
  this function initialises), orthonormalized, and their Ritz values written
  to mu[].  Returns the Lanczos dimension m; *converged reports the stop.
*/
static int
chiral_lanczos(int sigma, int nmodes, qpb_double eps, int mmax,
	       qpb_spinor_field *u, qpb_double *mu, int *converged)
{
  char sc = (sigma > 0) ? '+' : '-';
  qpb_spinor_field *lv = qpb_alloc(sizeof(qpb_spinor_field)*mmax);
  int n_alloc = 0;
  qpb_spinor_field av = qpb_spinor_field_init();
  qpb_spinor_field t1 = qpb_spinor_field_init();
  qpb_spinor_field t2 = qpb_spinor_field_init();
  qpb_double *alpha = qpb_alloc(sizeof(qpb_double)*mmax);
  qpb_double *beta  = qpb_alloc(sizeof(qpb_double)*mmax);

  /* v_0: random, projected onto the sector, normalized */
  lv[0] = qpb_spinor_field_init();
  n_alloc = 1;
  qpb_spinor_field_set_random(lv[0]);
  chiral_project(lv[0], lv[0], sigma, t1);
  qpb_double nrm;
  qpb_spinor_xdotx(&nrm, lv[0]);
  qpb_spinor_ax(lv[0], (qpb_complex){1./sqrt(nrm), 0.}, lv[0]);

  int m = 0;
  *converged = 0;
  qpb_double t = qpb_stop_watch(0);
  for(int j=0; j<mmax; j++)
    {
      apply_chiral(av, lv[j], sigma, t1, t2);
      if(j > 0)
	qpb_spinor_axpy(av, (qpb_complex){-beta[j-1], 0.}, lv[j-1], av);
      qpb_complex_double a;
      qpb_spinor_xdoty(&a, lv[j], av);
      alpha[j] = a.re;
      qpb_spinor_axpy(av, (qpb_complex){-alpha[j], 0.}, lv[j], av);

      /* full re-orthogonalization (twice), then remove any leakage out of
	 the sector */
      for(int pass=0; pass<2; pass++)
	for(int i=0; i<=j; i++)
	  {
	    qpb_complex_double c;
	    qpb_spinor_xdoty(&c, lv[i], av);
	    qpb_spinor_axpy(av, (qpb_complex){-c.re, -c.im}, lv[i], av);
	  }
      chiral_project(av, av, sigma, t1);

      qpb_double bb;
      qpb_spinor_xdotx(&bb, av);
      beta[j] = sqrt(bb);
      m = j+1;

      /* Ritz-residual test, every step up to m = 100, then every 10 */
      if(m >= nmodes+1 && (m <= 100 || m % 10 == 0))
	{
	  gsl_matrix *T = gsl_matrix_calloc(m, m);
	  for(int i=0; i<m; i++)
	    gsl_matrix_set(T, i, i, alpha[i]);
	  for(int i=0; i<m-1; i++)
	    {
	      gsl_matrix_set(T, i, i+1, beta[i]);
	      gsl_matrix_set(T, i+1, i, beta[i]);
	    }
	  gsl_vector *ev = gsl_vector_alloc(m);
	  gsl_matrix *S  = gsl_matrix_alloc(m, m);
	  gsl_eigen_symmv_workspace *ws = gsl_eigen_symmv_alloc(m);
	  gsl_eigen_symmv(T, ev, S, ws);
	  gsl_eigen_symmv_sort(ev, S, GSL_EIGEN_SORT_VAL_ASC);
	  gsl_eigen_symmv_free(ws);

	  qpb_double rmax = 0.;
	  for(int i=0; i<nmodes; i++)
	    {
	      qpb_double r = beta[j]*fabs(gsl_matrix_get(S, m-1, i));
	      if(r > rmax) rmax = r;
	    }
	  if(m % 20 == 0)
	    {
	      print("   sector %c: m = %4d, 1+mu_1 = %e, 1+mu_%d = %e, "
		    "max Ritz residual = %e (target = %e)\n",
		    sc, m, 1.+gsl_vector_get(ev, 0), nmodes,
		    1.+gsl_vector_get(ev, nmodes-1), rmax, eps);
	    }
	  gsl_matrix_free(T);
	  gsl_vector_free(ev);
	  gsl_matrix_free(S);
	  if(rmax < eps)
	    {
	      *converged = 1;
	      break;
	    }
	}

      /* next Lanczos vector */
      if(j < mmax-1)
	{
	  lv[j+1] = qpb_spinor_field_init();
	  n_alloc = j+2;
	  if(beta[j] < 1e-12)
	    {
	      /* invariant subspace: restart with a fresh random vector in the
		 sector, orthogonalized against the basis */
	      qpb_spinor_field_set_random(lv[j+1]);
	      chiral_project(lv[j+1], lv[j+1], sigma, t1);
	      for(int i=0; i<=j; i++)
		{
		  qpb_complex_double c;
		  qpb_spinor_xdoty(&c, lv[i], lv[j+1]);
		  qpb_spinor_axpy(lv[j+1], (qpb_complex){-c.re, -c.im},
				  lv[i], lv[j+1]);
		}
	      qpb_double nn;
	      qpb_spinor_xdotx(&nn, lv[j+1]);
	      qpb_spinor_ax(lv[j+1], (qpb_complex){1./sqrt(nn), 0.}, lv[j+1]);
	      beta[j] = 0.;
	    }
	  else
	    {
	      qpb_spinor_ax(lv[j+1], (qpb_complex){1./beta[j], 0.}, av);
	    }
	}
    }
  t = qpb_stop_watch(t);

  if(*converged)
    {
      print(" Chiral Lanczos, sector %c: converged, m = %d, t = %g sec\n",
	    sc, m, t);
    }
  else
    {
      error(" Chiral Lanczos, sector %c: WARNING -- hit max vectors = %d "
	    "without reaching eps = %e; the modes may be under-converged\n",
	    sc, mmax, eps);
      print(" Chiral Lanczos, sector %c: WARNING -- hit max vectors = %d "
	    "without reaching eps = %e; the modes may be under-converged\n",
	    sc, mmax, eps);
    }

  /* final eigendecomposition and lifting of the nmodes lowest Ritz vectors */
  gsl_matrix *T = gsl_matrix_calloc(m, m);
  for(int i=0; i<m; i++)
    gsl_matrix_set(T, i, i, alpha[i]);
  for(int i=0; i<m-1; i++)
    {
      gsl_matrix_set(T, i, i+1, beta[i]);
      gsl_matrix_set(T, i+1, i, beta[i]);
    }
  gsl_vector *ev = gsl_vector_alloc(m);
  gsl_matrix *S  = gsl_matrix_alloc(m, m);
  gsl_eigen_symmv_workspace *ws = gsl_eigen_symmv_alloc(m);
  gsl_eigen_symmv(T, ev, S, ws);
  gsl_eigen_symmv_sort(ev, S, GSL_EIGEN_SORT_VAL_ASC);
  gsl_eigen_symmv_free(ws);

  for(int c=0; c<nmodes; c++)
    {
      mu[c] = gsl_vector_get(ev, c);
      u[c] = qpb_spinor_field_init();
      qpb_spinor_field_set_zero(u[c]);
      for(int i=0; i<m; i++)
	qpb_spinor_axpy(u[c], (qpb_complex){gsl_matrix_get(S, i, c), 0.},
			lv[i], u[c]);
    }
  gsl_matrix_free(T);
  gsl_vector_free(ev);
  gsl_matrix_free(S);

  /* modified Gram-Schmidt: enforce U^dag U = 1 within the sector */
  for(int a=0; a<nmodes; a++)
    {
      for(int b=0; b<a; b++)
	{
	  qpb_complex_double c;
	  qpb_spinor_xdoty(&c, u[b], u[a]);
	  qpb_spinor_axpy(u[a], (qpb_complex){-c.re, -c.im}, u[b], u[a]);
	}
      qpb_double n2;
      qpb_spinor_xdotx(&n2, u[a]);
      qpb_spinor_ax(u[a], (qpb_complex){1./sqrt(n2), 0.}, u[a]);
    }

  for(int i=0; i<n_alloc; i++)
    qpb_spinor_field_finalize(lv[i]);
  free(lv);
  qpb_spinor_field_finalize(av);
  qpb_spinor_field_finalize(t1);
  qpb_spinor_field_finalize(t2);
  free(alpha);
  free(beta);
  return m;
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

  /* am values for the mass checks: one chiral build serves every mass */
  qpb_double masses[CHIRAL_MAX_MASSES];
  int n_masses = 0;
  {
    char masses_str[QPB_MAX_STRING];
    strncpy(masses_str, qpb_parse("Masses"), QPB_MAX_STRING-1);
    masses_str[QPB_MAX_STRING-1] = '\0';
    for(char *tok = strtok(masses_str, " \t"); tok != NULL;
	tok = strtok(NULL, " \t"))
      {
	if(n_masses == CHIRAL_MAX_MASSES)
	  {
	    error("Masses: at most %d values, quiting\n", CHIRAL_MAX_MASSES);
	    exit(QPB_PARAMETERS_ERROR);
	  }
	if(sscanf(tok, "%lf", &masses[n_masses])!=1)
	  {
	    error("error parsing for %s\n", "Masses");
	    exit(QPB_PARSER_ERROR);
	  }
	n_masses++;
      }
    if(n_masses == 0)
      {
	error("error parsing for %s\n", "Masses");
	exit(QPB_PARSER_ERROR);
      }
  }

  qpb_double c_sw;
  if(sscanf(qpb_parse("c_sw"), "%lf", &c_sw)!=1)
    {
      error("error parsing for %s\n",
	    "c_sw");
      exit(QPB_PARSER_ERROR);
    }

  /* Inner multi-shift CG (sign function) parameters. NOTE: the
     stopping criterion is on the *squared* relative residual. Since
     Lanczos relies on H being Hermitian, a tight value is advisable */
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

  /* Chiral Lanczos on A± = P± g5 s P±, one sector at a time */
  int chiral_modes;
  if(sscanf(qpb_parse("Chiral modes"), "%d", &chiral_modes)!=1)
    {
      error("error parsing for %s\n",
	    "Chiral modes");
      exit(QPB_PARSER_ERROR);
    }
  if(chiral_modes < 1 || chiral_modes > CHIRAL_MAX_MODES)
    {
      error("Chiral modes = %d out of range [1, %d], quiting\n",
	    chiral_modes, CHIRAL_MAX_MODES);
      exit(QPB_PARAMETERS_ERROR);
    }

  qpb_double chiral_epsilon;
  if(sscanf(qpb_parse("Chiral Lanczos epsilon"), "%lf", &chiral_epsilon)!=1)
    {
      error("error parsing for %s\n",
	    "Chiral Lanczos epsilon");
      exit(QPB_PARSER_ERROR);
    }
  if(chiral_epsilon <= 0.)
    {
      error("Chiral Lanczos epsilon = %g must be positive, quiting\n",
	    chiral_epsilon);
      exit(QPB_PARAMETERS_ERROR);
    }

  int chiral_max_vectors;
  if(sscanf(qpb_parse("Chiral Lanczos max vectors"), "%d",
	    &chiral_max_vectors)!=1)
    {
      error("error parsing for %s\n",
	    "Chiral Lanczos max vectors");
      exit(QPB_PARSER_ERROR);
    }
  if(chiral_max_vectors < chiral_modes+1
     || chiral_max_vectors > QPB_KL_DEFL_MAX_VECTORS)
    {
      error("Chiral Lanczos max vectors = %d out of range [%d, %d], quiting\n",
	    chiral_max_vectors, chiral_modes+1, QPB_KL_DEFL_MAX_VECTORS);
      exit(QPB_PARAMETERS_ERROR);
    }

  /* Sign-function deflation build (inside qpb_overlap_kl_pfrac_init): the
     Lanczos on X^2 whose low modes are treated exactly */
  qpb_double Lanczos_epsilon;
  if(sscanf(qpb_parse("Lanczos epsilon"), "%lf", &Lanczos_epsilon)!=1)
    {
      error("error parsing for %s\n",
	    "Lanczos epsilon");
      exit(QPB_PARSER_ERROR);
    }

  int Lanczos_max_vectors;
  if(sscanf(qpb_parse("Lanczos max vectors"), "%d", &Lanczos_max_vectors)!=1)
    {
      error("error parsing for %s\n",
	    "Lanczos max vectors");
      exit(QPB_PARSER_ERROR);
    }
  if(Lanczos_max_vectors < 2 || Lanczos_max_vectors > QPB_KL_DEFL_MAX_VECTORS)
    {
      error("Lanczos max vectors = %d out of range [2, %d], quiting\n",
	    Lanczos_max_vectors, QPB_KL_DEFL_MAX_VECTORS);
      exit(QPB_PARAMETERS_ERROR);
    }

  /* Optional: "Deflation = off" skips the low-mode deflation build entirely
     (whatever QPB_DEFL_K is); "on" or no key keeps the default behaviour. */
  int defl_off = 0;
  char *defl_switch = qpb_parse_optional("Deflation");
  if(defl_switch != NULL)
    {
      char defl_word[QPB_MAX_STRING];
      if(sscanf(defl_switch, "%s", defl_word)!=1)
	{
	  error("error parsing for %s\n", "Deflation");
	  exit(QPB_PARSER_ERROR);
	}
      if(strcmp(defl_word, "off") == 0)
	defl_off = 1;
      else if(strcmp(defl_word, "on") != 0)
	{
	  error("%s: option should be one of: ", "Deflation");
	  error("%s, ", "on");
	  error("%s\n", "off");
	  exit(QPB_PARSER_ERROR);
	}
    }
  if(defl_off)
    qpb_overlap_kl_pfrac_disable_deflation();

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
  print(" Masses =");
  for(int i=0; i<n_masses; i++)
    {
      print(" %g", masses[i]);
    }
  print("\n");
  print(" Clover param = %g\n", c_sw);
  print(" BC in time = %g\n", timebc);
  print(" Solver epsilon = %e\n", epsilon);
  print(" Max solver iters = %d\n", max_iters);
  print(" Chiral modes = %d\n", chiral_modes);
  print(" Chiral Lanczos epsilon = %e\n", chiral_epsilon);
  print(" Chiral Lanczos max vectors = %d\n", chiral_max_vectors);
  print(" Lanczos epsilon = %e\n", Lanczos_epsilon);
  print(" Lanczos max vectors = %d\n", Lanczos_max_vectors);
  if(defl_off)
    {
      print(" Deflation = off\n");
    }
  switch(which_dslash_op)
    {
    case QPB_DSLASH_BRILLOUIN:
      print(" Dslash operator is Brillouin\n");
      break;
    case QPB_DSLASH_STANDARD:
      print(" Dslash operator is Standard\n");
      break;
    }
  print(" KL class = %d\n", 11);
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

  /* Only the sign function is applied, so the overlap mass passed here is
     irrelevant (0) and the preconditioner is switched off (prec_max_iter =
     0); its tolerances are inert.  D at each am is assembled below. */
  qpb_overlap_kl_pfrac_init(solver_arg_links, clover_term, kl_class, kl_iters,
                          rho, c_sw, 0., scaling_factor,
                          epsilon, epsilon, max_iters,
                          epsilon, 0,
                          Lanczos_epsilon, Lanczos_max_vectors);

  /* ---- 1. Hermiticity of s:  |<x, s y> - <s x, y>| / (|x| |y|) ---- */
  {
    qpb_spinor_field x  = qpb_spinor_field_init();
    qpb_spinor_field y  = qpb_spinor_field_init();
    qpb_spinor_field sx = qpb_spinor_field_init();
    qpb_spinor_field sy = qpb_spinor_field_init();
    qpb_spinor_field_set_gaussian(x);
    qpb_spinor_field_set_gaussian(y);
    apply_sign(sx, x);
    apply_sign(sy, y);
    qpb_complex_double c1, c2;
    qpb_double nx, ny;
    qpb_spinor_xdoty(&c1, x, sy);
    qpb_spinor_xdoty(&c2, sx, y);
    qpb_spinor_xdotx(&nx, x);
    qpb_spinor_xdotx(&ny, y);
    print(" Hermiticity check: |<x, s y> - <s x, y>| / (|x| |y|) = %e\n",
	  sqrt((c1.re-c2.re)*(c1.re-c2.re) + (c1.im-c2.im)*(c1.im-c2.im))
	  / sqrt(nx*ny));
    qpb_spinor_field_finalize(x);
    qpb_spinor_field_finalize(y);
    qpb_spinor_field_finalize(sx);
    qpb_spinor_field_finalize(sy);
  }

  /* ---- 2. chiral Lanczos, one sector at a time ---- */
  qpb_spinor_field *u = qpb_alloc(sizeof(qpb_spinor_field)*chiral_modes);
  qpb_spinor_field *G = qpb_alloc(sizeof(qpb_spinor_field)*chiral_modes);
  qpb_double *mu = qpb_alloc(sizeof(qpb_double)*chiral_modes);
  qpb_spinor_field t1  = qpb_spinor_field_init();
  qpb_spinor_field t2  = qpb_spinor_field_init();
  qpb_spinor_field S1  = qpb_spinor_field_init();
  qpb_spinor_field S2  = qpb_spinor_field_init();
  qpb_spinor_field Du  = qpb_spinor_field_init();
  qpb_spinor_field DDu = qpb_spinor_field_init();

  for(int sec=0; sec<2; sec++)
    {
      int sigma = (sec == 0) ? +1 : -1;
      char sc = (sigma > 0) ? '+' : '-';
      int converged;
      print("\n Chiral Lanczos, sector %c...\n", sc);
      chiral_lanczos(sigma, chiral_modes, chiral_epsilon, chiral_max_vectors,
		     u, mu, &converged);

      /* explicit check of each mode, with G = g5 s u kept for the mass
	 checks:  A u = P g5 s u = P G  (u lies in the sector) */
      for(int c=0; c<chiral_modes; c++)
	{
	  G[c] = qpb_spinor_field_init();
	  qpb_gamma5_sign_function_of_X_pfrac(G[c], u[c]);
	  chiral_project(t1, G[c], sigma, t2);              /* t1 = A u   */
	  qpb_complex_double z;
	  qpb_spinor_xdoty(&z, u[c], t1);                   /* mu = <u, A u> */
	  qpb_spinor_axpy(t1, (qpb_complex){-z.re, 0.}, u[c], t1);
	  qpb_double rr;
	  qpb_spinor_xdotx(&rr, t1);
	  print(" Chiral sector %c, mode %2d: 1+mu = %+.12e, "
		"||A u - mu u|| = %e\n", sc, c, 1.+z.re, sqrt(rr));
	}

      /* ---- 3. mass checks:  D = a + b g5 s,  D^dag = a + b s g5 ----
	 With G = g5 s u, S1 = s g5 u and S2 = s g5 G (two sign applications
	 per mode), every am is plain vector arithmetic:
	   D u       = a u + b G
	   D^dag D u = a D u + b (a S1 + b S2)                                */
      for(int c=0; c<chiral_modes; c++)
	{
	  qpb_spinor_gamma5(t1, u[c]);
	  apply_sign(S1, t1);                              /* S1 = s g5 u */
	  qpb_spinor_gamma5(t1, G[c]);
	  apply_sign(S2, t1);                              /* S2 = s g5 G */
	  for(int im=0; im<n_masses; im++)
	    {
	      qpb_double a = rho + 0.5*masses[im];
	      qpb_double b = rho - 0.5*masses[im];
	      qpb_spinor_axpby(Du, (qpb_complex){a, 0.}, u[c],
			       (qpb_complex){b, 0.}, G[c]);       /* D u */
	      qpb_spinor_axpby(DDu, (qpb_complex){a*b, 0.}, S1,
			       (qpb_complex){b*b, 0.}, S2);
	      qpb_spinor_axpy(DDu, (qpb_complex){a, 0.}, Du, DDu); /* D^dag D u */
	      qpb_double lam2;
	      qpb_spinor_xdotx(&lam2, Du);                  /* <u, D^dag D u> */
	      qpb_spinor_axpy(DDu, (qpb_complex){-lam2, 0.}, u[c], DDu);
	      qpb_double rr;
	      qpb_spinor_xdotx(&rr, DDu);
	      print(" Mass check, sector %c, mode %2d, am = %g: lambda^2 = %.12e, "
		    "||D^dag D u - lambda^2 u|| = %e\n",
		    sc, c, masses[im], lam2, sqrt(rr));
	    }
	}

      for(int c=0; c<chiral_modes; c++)
	{
	  qpb_spinor_field_finalize(u[c]);
	  qpb_spinor_field_finalize(G[c]);
	}
    }
  free(u);
  free(G);
  free(mu);
  qpb_spinor_field_finalize(t1);
  qpb_spinor_field_finalize(t2);
  qpb_spinor_field_finalize(S1);
  qpb_spinor_field_finalize(S2);
  qpb_spinor_field_finalize(Du);
  qpb_spinor_field_finalize(DDu);

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
