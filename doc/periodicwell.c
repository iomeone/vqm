/*
  periodicwell.c
*/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define ABS(a)     (((a) < 0) ? -(a) : (a))
/* equivalent of Fortran nint intrinsic */
#define NINT(a) ((a) >= 0.0 ? (int)((a)+0.5) : (int)((a)-0.5))

/* void cft_( double*, double*, int*, int*, int*, int*);
   void dsyev_( char*, char*, int*, double*, int*, double; int*, int*); */

void main()
{
  /* 
     Band structure of a 1d model periodic system (Kronig-Penney)
     Expansion on a plane-wave basis set and diagonalization
     Units: hbar^2/2m = 1
     Requires lapack dsyev, cern library cft
  */
  
  static double pi = 3.14159265358979;
  int n, npw;
  double v0, a, b, ecut, k, x;
  double *g, *e, *h, *wrk;
  double *vr, *vi;
  int i, j, ij, m, lwork, info;
  int nfft, ifft, signfft=-1;
  char *V = "V";
  char *U = "U";
  FILE *out;
  
  /*  Input data
      Potential well: V(x)=-V_0 for |x|<b/2, V(x)=0 for |x|>b/2
      Periodicity:    V(x+a)=V(x)
  */
  fprintf(stdout, "Parameters for potential well: V_0, a, b > ");
  scanf("%lf %lf %lf",&v0, &a, &b);
  if ( v0 <= 0 ||  a <= 0 || b <= 0 || a <= b )  { 
    fprintf(stderr, "wrong input parameters\n");
    exit(1);
  }
  fprintf (stdout,"   V_0=%f, a0=%f, b=%f\n",v0,a,b);
  /*
    Plane waves basis set: G_n=n*2pi/a, \hbar^2/2m*G^2 < Ecut
  */
  fprintf(stdout, "Cutoff for plane waves: ecut > ");
  scanf("%lf",&ecut);
  if ( ecut <= 0) {
    fprintf(stderr, "wrong input parameter\n");
    exit(1);
  }
  /*
    Number of plane waves
  */
  npw = (int) ( sqrt ( ecut/pow( 2.0*pi/a, 2) ) + 0.5 );
  npw = 2*npw+1;
  fprintf (stdout,"   ecut = %f,  # of PWs=%d\n",ecut, npw);
  /*
    Assign values of  G_n: n=0,+1,-1,+2,-2, etc
  */
  g  = (double *) malloc ( npw * sizeof (double) );
  e  = (double *) malloc ( npw * sizeof (double) );
  h  = (double *) malloc (npw*npw*sizeof (double) );
  wrk= (double *) malloc (3*npw* sizeof (double) );

  g[0] = 0.0;
  for ( i = 1; i < npw; i+=2 ) {
    g[i  ] = (i+1)*pi/a;
    g[i+1] =-(i+1)*pi/a;
  }
  for ( ij = 0; ij < npw*npw; ++ij ) {
    h[ij] = 0.0;
  }
  /*    Compute V(G) with FFT
        nfft = number of points in real and reciprocal space
               must be at least four times the number of plane waves
               because we need V(G) for G=G_i-G_j
  */
  nfft = 4*npw;
  fprintf (stdout,"FFT dimension: min=%d   Your value > ",nfft);
  scanf("%d",&nfft);
  vr = ( double *) malloc ( nfft* sizeof (double) );
  vi = ( double *) malloc ( nfft* sizeof (double) );
  /*
     compute v(x_i), x_i=(i-1)*dx, dx=a/nfft
  */
  for ( i = 0; i < nfft; i++ ) {
    x =i*(a/nfft);
    if ( x <= b/2.0 ||  x > (a-b/2.0) ) {
      vr[i] = -v0; }
    else {
      vr[i] = 0.0;
      /* 
	 imaginary part is zero in this case, for both real and reciprocal space
      */
      vi[i] = 0.0;
    }
  }

  cft_ ( vr, vi, &nfft, &nfft, &nfft, &signfft ) ;
    for ( i = 0; i < nfft; i++ ) { vr[i] = vr[i]/nfft ; }
  /*
    Loop on k-vectors: k runs from -pi/a to pi/a
  */
  out = fopen("bands.out", "w");
  n = 20;
  for ( m =-n; m <= n; m++ ) {
    k = m*pi/n/a;
    /*
      Assign values of the matrix elements of the hamiltonian 
      on the plane wave basis
    */
    ij = 0;
    for ( i = 0; i < npw; ++i ) {
      for ( j = 0; j < npw; ++j ) {
      /* NOTA BENE: the matrix h is a vector in the calling program,
         while dsyev expects a (pointer to a) fortran matrix.
         A fortran matrix is "simulated" in the following way:
         if h[ij] is the C array and  h(i,j) is a N*M Fortran matrix,
         h[ij] == h(i,j), where  ij = (j-1)*N + (i-1)  (i=1,N, j=1,M)
      */

	if ( i == j ) {
	  h[ij++] = (k+g[i])*(k+g[i]) + vr[0];
	  /* vr[0] = v(G=0) = -v0/a*b */
	} else {
	  /* ifft points to the component G=g(j)-g(i) */
	  ifft = NINT ( (g[j]-g[i])*a/pi ) / 2 ;
	  if ( ifft < 0 ) { ifft += nfft ; }
	  h[ij++] = vr[ifft];
	  /* h(i,j) = -v0/a * sin((g(j)-g(i))*b/2) / (g(j)-g(i))*2 */
	}
      }
    }
    /*
      Solution [expansion coefficients are stored into h(i,j)
      j=basis function index, i= eigenvalue index]
      (beware fortran-C reversed index convention!)
    */
    lwork = 3*npw;
    dsyev_ ( V, U, &npw, h, &npw, e, wrk, &lwork, &info ); 
    if ( info != 0) {
      fprintf(stderr, "H-matrix diagonalization failed\n");
      exit(1);
    }
    printf("k=%f, Lowest eigenvalues: %f %f %f\n",k,e[0],e[1],e[2]);
    /*
      Write to output file the band dispersion e(k)
    */
    fprintf(out,"%12.6f %10.6f %10.6f %10.6f\n",k,e[0],e[1],e[2]);
  }
  fclose(out);
}
