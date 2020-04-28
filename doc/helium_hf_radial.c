/* helium_hf_radial.c */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define ABS(a)	   (((a) < 0) ? -(a) : (a))
#define MIN(a,b)   (((a) < (b)) ? (a) : (b))
#define MAX(a,b)   (((a) > (b)) ? (a) : (b))

static double pi=3.14159265358979;

main()
{
  /* Subroutines */
    extern int do_mesh(int, double, double, double, double, 
		       double *, double *, double *);
    extern void init_pot (double, int, double *, double *);
    extern double solve_sheq(int, int, double, int, double,
			     double *, double *, double *, double *, 
			     double *);
    extern void rho_of_r  (int, double *, double*, double*, double *);
    extern double v_of_rho (int, double, double *, double *,
			    double *, double *);
  /* Variables */
    double *r, *r2, *sqr;
    double *y, *vpot, *rho, *vscf, *vhx, *deltav;
    int mesh;
    double zeta, rmax, xmin, dx, zmesh, e, de, beta, tol;
    double ehx, evion, ekin, etot, etot1;
    int i, l, n, iter;
    FILE *out;

/* read atomic charge and other parameters */
    fprintf(stdout,"Hartree-Fock calculation for Helium-like atoms\n");
    fprintf(stdout, " ATOMIC CHARGE = ");
    scanf("%lf",&zeta);
    fprintf(stdout, " zeta = %f\n", zeta);
    if ( zeta < 2.0) {
      fprintf(stderr, "zeta should be > 2\n");
      exit (1);
    }
    fprintf(stdout, " MIXING PARAMETER beta [0.0-1.0] = ");
    scanf("%lf",&beta);
    fprintf(stdout, " beta = %f\n", beta);
    if ( beta <= 0.0 || beta > 1.0) {
      fprintf(stderr, "beta out of range\n");
      exit (1);
    }
    fprintf(stdout, " SCF ACCURACY (Ry) = ");
    scanf("%lf",&tol);
    fprintf(stdout, " tol= %f\n", tol);
    if ( tol <= 0.0 ) {
      fprintf(stderr, "tol should be strictly positive\n");
      exit (1);
    }

    /* initialize logarithmic mesh */

    zmesh = zeta;
    rmax = 100.;
    xmin = -6.0;
    dx = 0.01 ;
/* number of grid points */
    mesh = (int ) ((log(zmesh * rmax) - xmin) / dx);

    r = malloc( (mesh+1) * sizeof(double) );
    r2= malloc( (mesh+1) * sizeof(double) );
    sqr=malloc( (mesh+1) * sizeof(double) );

    do_mesh ( mesh, zmesh, xmin, dx, rmax, r, sqr, r2);

    /* initialize the potential */

    vpot = malloc( (mesh+1) * sizeof(double) );
    init_pot(zeta, mesh, r, vpot);

    /*  The GS configuration for helium-like atoms is (1s)**2  */

    n = 1;
    l = 0;

    vscf = malloc( (mesh+1) * sizeof(double) );
    for (i = 0; i <= mesh; i++) {
      vscf[i] = vpot[i];
    }

    /* SCF cycle */
    y   = malloc( (mesh+1) * sizeof(double) );
    rho = malloc( (mesh+1) * sizeof(double) );
    vhx = malloc( (mesh+1) * sizeof(double) );
    deltav = malloc( (mesh+1) * sizeof(double) );

    /* the variational correction "de" is used to check for self-consistency */
    de = 2.*tol;

    for (iter = 1; ABS(de) > tol; ++iter){
      fprintf(stdout, "####################################################\n");
      fprintf(stdout, " SCF iteration # %d\n",iter);

    /* solve schroedinger equation in radial coordinates by Numerov method */

    e = solve_sheq(n, l, zeta, mesh, dx, r, sqr, r2, vscf, y);

    /* calculate the charge density from the wfc */

    rho_of_r ( mesh, r, r2, y, rho) ;

    /* calculate the Hartree + Exchange potential and energy */

    ehx = v_of_rho ( mesh, dx, r, r2, rho, vhx );
 
    /* calculate the kinetic energy and the energy in the external potential */

    evion = 0.0;
    ekin = 2.0 * e;
    for (i = 0; i <= mesh; i++){
      evion = evion + vpot[i] * rho[i] * 4.*pi * r2[i] * r[i] * dx;
      ekin  = ekin  - vscf[i] * rho[i] * 4.*pi * r2[i] * r[i] * dx;
    }
    de = 0.0;

    for (i = 0; i <= mesh; i++){
      deltav[i] = vpot[i] + vhx[i] - vscf[i];
      vscf[i] = vscf[i] + beta * deltav[i];
      de = de + deltav[i] * rho[i] * 4.*pi * r2[i] * r[i] * dx;
    }

    /* write out the eigenvalue energy 
       to be compared with the external potential */

    etot = 2.0 * e - ehx + de;
    etot1= ekin + evion + ehx;
    fprintf (stdout, "eigenvalue  = %lf\n" , e);
    fprintf (stdout, "Eigenvalue energy      %15.6f\n" , 2*e);
    fprintf (stdout, "Kinetic energy         %15.6f\n", ekin);
    fprintf (stdout, "External pot. energy   %15.6f\n", evion);
    fprintf (stdout, "Hartree+Exch. energy   %15.6f\n", ehx);
    fprintf (stdout, "Variational correct.   %15.6f\n",de);
    fprintf (stdout, "Total energy           %15.6f %15.6f\n", etot, etot1);
    fprintf (stdout, "Virial check           %15.6f\n",-(evion+ehx)/ekin);
    }

    fprintf (stdout, " SCF Convergence has been achieved\n");
    fprintf (stdout, " compute additional single-particle states\n");

    /* write to file pot.out */

    out = fopen("pot.out", "w");
    for (i = 0; i <= mesh; i++ ){
      fprintf(out, "%lf %lf %lf %lf\n",r[i],vpot[i],vhx[i],vscf[i]);
    }
    fclose(out);

      /* open wfc file */

    out = fopen("wfc.out", "w");

 L1:
      /* read principal and angular quantum numbers */

    fprintf(stdout, " n,l >> " );
    scanf("%d %d",&n, &l);
    if (n < 1) {
      fclose(out);
      free(deltav);
      free(vhx);
      free(rho);
      free(y); 
      free(vscf); 
      free(vpot);
      free(sqr);
      free(r2);
      free(r);
      exit (0);
    }
    if (n < l + 1) {
      fprintf(stderr, "error in main: n.lt.l -> wrong number of nodes\n");
      exit (1);
    }

/* solve the schroedinger equation in radial coordinates by Numerov method */

    e = solve_sheq(n, l, zeta, mesh, dx, r, sqr, r2, vscf, y);

    fprintf (stdout, "eigenvalue %d %d = %16.8e\n", n, l, e );

    for (i = 0; i <= mesh; ++i) {
      fprintf(out, "%14.6e %14.6e %14.6e %14.6e\n", 
	      r[i], y[i] / sqr[i], y[i] * sqr[i], e);
    }
    fprintf(out, "\n\n"); 
    goto L1;
} /* MAIN__ */

/* --------------------------------------------------------------------- */
/* Subroutine */ double solve_sheq(int n, int l, double zeta, int mesh, 
				   double dx, double *r, double *sqr, 
				   double *r2, double *vpot, double *y ) 

{
    /* Local variables */
     int i, j;
     double e, de, fac;
     int icl, kkk;
     double x2l2, elw, eup, ddx12, norm;
     int nodes;
     double sqlhf, ycusp, dfcusp;
     int ncross;
     double *f;

/* --------------------------------------------------------------------- */

/* solve the schroedinger equation in radial coordinates on a 
   logarithmic grid by Numerov method - atomic (Ry) units */


    ddx12 = dx * dx / 12.;
/* Computing 2nd power */
    sqlhf = (l + 0.5) * (l + 0.5);
    x2l2 = (double) (2*l+ 2);

/* set initial lower and upper bounds to the eigenvalue */

    eup = vpot[mesh];
    elw = eup;
    for (i = 0; i <= mesh; ++i) {
/*      if ( elw > sqlhf / r2[i] + vpot[i] ) 
	elw = sqlhf / r2[i] + vpot[i] ; */
      elw = MIN ( elw, sqlhf / r2[i] + vpot[i] );
    }
    if (eup - elw < 1e-10) {
      fprintf (stderr, "%25.16e 25.16e\n", eup, elw);
      fprintf (stderr, "solve_sheq: lower and upper bounds are equal\n");
      exit(1);
    }
    e = (elw + eup) * .5;
    f = malloc( (mesh+1) * sizeof(double) );
    kkk = 0;
L1:
/* this is the entry point for the solution at fixed energy */
    ++kkk;

/* set up the f-function and determine the position of its last */
/* change of sign */
/* f < 0 (approximately) means classically allowed   region */
/* f > 0         "         "        "      forbidden   " */

    icl = -1;
    f[0] = ddx12 * (sqlhf + r2[0] * (vpot[0] - e));
    for (i = 1; i <= mesh; ++i) {
	f[i] = ddx12 * (sqlhf + r2[i] * (vpot[i] - e));
/* beware: if f(i) is exactly zero the change of sign is not observed */
/* the following line is a trick to prevent missing a change of sign */
/* in this unlikely but not impossible case: */
	if (f[i] == 0.) {
	    f[i] = 1e-20;
	}
	if (f[i] != copysign(f[i], f[i - 1])) {
	    icl = i;
	}
    }
    if (icl < 0) {
      fprintf (stderr, "solve_sheq: no classical turning point");
      exit(1);
    }
    if (icl < 0 || icl >= mesh) {
/*
      fprintf (stderr, "%4d %4d\n", icl, mesh);
      fprintf (stderr, "error in solve_sheq: last change of sign too far");
      exit(1); 
*/
      eup = e;
      e = (eup + elw) * .5;
      goto L1;
    }

/* f function as required by numerov method */

    for (i = 0; i <= mesh; ++i) {
	f[i] = 1. - f[i];
	y[i] = 0.;
    }

/* determination of the wave-function in the first two points */

    nodes = n - l - 1;
    y[0] = pow (r[0], l+1) * (1. - zeta * 2. * r[0] / x2l2) / sqr[0];
    y[1] = pow (r[1], l+1) * (1. - zeta * 2. * r[1] / x2l2) / sqr[1];

/* outward integration, count number of crossings */

    ncross = 0;
    for (i = 1; i <= icl-1; ++i) {
	y[i + 1] = ((12. - f[i] * 10.) * y[i] - f[i - 1] * y[i - 1])
		 / f[i + 1];
	if (y[i] != copysign(y[i],y[i+1]) ) {
	    ++ncross;
	}
    }
    fac = y[icl];

/* check number of crossings */

    if (ncross != nodes) {
	if (kkk > 100) {
	  fprintf(stderr, "%4d %4d %4d %4d %16.8e %16.8e %16.8e\n", 
		  kkk, ncross, nodes, icl,  e, elw, eup);
	  fprintf(stderr, " error in solve_sheq: too many iterations\n");
	  exit (1);
	}
	if (ncross > nodes) {
	    eup = e;
	} else {
	    elw = e;
	}
	e = (eup + elw) * .5;
	goto L1;
    }

/* determination of the wave-function in the last two points */
/* assuming y(mesh+1) = 0 and y(mesh) = dx */

    y[mesh] = dx;
    y[mesh - 1] = (12. - f[mesh] * 10.) * y[mesh] / f[mesh - 1];

/* inward integration */

    for (i = mesh - 1; i >= icl+1; --i) {
	y[i - 1] = ((12. - f[i] * 10.) * y[i] - f[i + 1] * y[i + 1])
		 / f[i - 1];
	if (y[i - 1] > 1e10) {
	    for (j = mesh; j >= i-1; --j) {
		y[j] /= y[i - 1];
	    }
	}
    }

/* rescale function to match at the classical turning point (icl) */

    fac /= y[icl];
    for (i = icl; i <= mesh; ++i) {
	y[i] *= fac;
    }

/* normalize on the segment */

    norm = 0.;
    for (i = 1; i <= mesh; ++i) {
	norm += y[i] * y[i] * r2[i] * dx;
    }
    norm = sqrt(norm);
    for (i = 0; i <= mesh; ++i) {
	y[i] /= norm;
    }

/* find the value of the cusp at the matching point (icl) */

    i = icl;
    ycusp = (y[i - 1] * f[i - 1] + f[i + 1] * y[i + 1] + f[i] * 10. 
	    * y[i]) / 12.;
    dfcusp = f[i] * (y[i] / ycusp - 1.);

/* eigenvalue update using perturbation theory */

    de = dfcusp / ddx12 * ycusp * ycusp * dx;
    if (de > 0.) {
	elw = e;
    }
    if (de < 0.) {
	eup = e;
    }

/* prevent e to go out of bounds, i.e. e > eup or e < elw */
/* (might happen far from convergence) */

    e = e + de;
    e = MIN (e,eup);
    e = MAX (e,elw);
    /* if ( e > eup ) e=eup;
       if ( e < elw ) e=elw; */
    if (kkk > 100) {
      fprintf(stderr, "%4d %16.8e %16.8e\n", kkk, e, de);
      fprintf(stderr, " error in solve_sheq: too many iterations\n");
      exit (1);
    }
    if (ABS(de) > 1e-10) {
	goto L1;
    }
/* ---- convergence has been achieved ----- */
    fprintf(stdout, "convergence achieved at iter # %4d, de = %16.8e\n",
	    kkk, de);
    free(f);
    return e;
} /* solve_sheq__ */

/* -------------------------------------------------------------------- */
/* Subroutine */ int do_mesh (int mesh, double zmesh, double xmin, 
	double dx, double rmax, 
	double *r, double *sqr, double *r2)
{

    /* Builtin functions */
    double log(double), exp(double), sqrt(double);

    /* Local variables */
    int i;
    double x;

/* -------------------------------------------------------------------- */

/* initialize grid */

    for (i = 0; i <= mesh; ++i ) {
	x = xmin + dx * i;
	r[i] = exp(x) / zmesh;
	sqr[i] = sqrt(r[i]);
	r2[i] = r[i] * r[i];
    }
    fprintf(stdout, " radial mesh information:\n");
    fprintf(stdout, " dx   = %12.6f", dx);
    fprintf(stdout, ", xmin = %12.6f", xmin);
    fprintf(stdout, ", zmesh =%12.6f\n", zmesh);
    fprintf(stdout, " mesh = %5d", mesh);
    fprintf(stdout, ", r(0) = %12.6f",  r[0]);
    fprintf(stdout, ", r(mesh) = %12.6f\n", r[mesh]);
    return mesh; 
} /* do_mesh */

/* -------------------------------------------------------------------- */
/* Subroutine */ void init_pot(double zeta, int mesh, double *r, 
			      double *vpot)
{
    /* Local variables */
    static int i;
    FILE *out;

/* -------------------------------------------------------------------- */

/* initialize potential */

    out = fopen("pot.out","w");
    for (i = 0; i <= mesh; ++i) {
	vpot[i] = -2 * zeta / r[i];
	fprintf(out, "%16.8e %16.8e\n", r[i], vpot[i]);
    }
    fclose(out);
    return;
} /* init_pot */

/* -------------------------------------------------------------------- */
/* Subroutine */ void rho_of_r (int mesh, double *r, double *r2,
				double *y, double *rho)
{
  /*  compute the charge density of a He-like atom  */
  
  double fpi, nelec=2.0;
  int i;
  
  fpi = 4.0*pi;
  for (i = 0; i <= mesh ; i++) {
    rho[i] = nelec * pow ( y[i], 2) * r[i] / (fpi*r2[i]);
  }
} /* rho_of_r */

/* -------------------------------------------------------------------- */
/* Subroutine */ double v_of_rho (int mesh, double dx, double *r,
				  double *r2, double *rho, double *vhx)
{
  /* compute the Hartree + Exchange potential for an He-like atom
     Exchange cancels exactly half of the Hartree result for both 
     potential and energy
  */

  static double e2=2.0;
  double fpi;
  int i;
  double charge, ehx;

  /* calculate the Hartree potential and energy by integrating the 
     electric field generated by the electronic charge. 
     This is done in 2 steps
     
     1) calculate the charge inside a sphere and fill vhx with the 
     electric field generated by this charge
  */
  charge = 0.0;
  fpi = 4.0*pi;
  for (i = 0; i <= mesh ; i++) {
    charge = charge + rho[i] * fpi * r2[i] * r[i] * dx;
    vhx[i] = e2*charge/r2[i];
  }

  /*  (the total charge is written in output as a check) */
  fprintf(stdout, "Total charge = %lf\n", charge);
  /*
    2) integrate the electric field from +\infty to r to get the potential
    and integrate V_{Hartree}*rho to get the energy
  */
  ehx = 0.0;
  vhx[mesh] = e2 * charge / r[mesh];
  for (i = mesh-1; i >= 0 ; i--) {
    vhx[i] = vhx[i+1] + vhx[i] * r[i] * dx;
    ehx = ehx + vhx[i] * rho[i] * fpi * r2[i] * r[i] * dx;
  }
  ehx = ehx/2.0; 
  /*
    Exchange cancels exactly half of the Hartree result for both 
    potential and energy
  */
  ehx = 0.5 * ehx; 
  for (i = 0; i <= mesh ; i++) {
    vhx[i] = 0.5 * vhx[i];
  }
  return ehx;
} /* v_of_rho */
