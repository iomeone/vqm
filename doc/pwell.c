/*
  pwell.c
*/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

main()
{
  /* 
     lowest energy levels of the finite (symmetric) potential well
     via expansion on a plane-wave basis set and diagonalization
     Units: hbar^2/2m = 1
     Requires lapack dsyev
  */
  
  static double pi = 3.14159265358979;
  int n, npw;
  double v0, a, b;
  double x, dx, norm, prob, fr, fi;
  double *h, *kn, *e, *wrk;
  int i, j, ij, nr, lwork, info;
  char *V = "V";
  char *U = "U";
  FILE *out;
  
  /*  Input data
      Potential well: V(x)=-V_0 for |x|<b/2, V(x)=0 for |x|>b/2
  */
  fprintf(stdout, "Parameters for potential well: V_0, b > ");
  scanf("%lf %lf",&v0, &b);
  if ( v0 <= 0 ||  b <= 0)  { 
    fprintf(stderr, "wrong input parameters\n");
    exit(1);
  }
  fprintf (stdout,"   V_0, b = %f %f\n",v0,b);
  /*
    Plane waves between -a/2<x<a/2, k_i=+-2*pi*i/a, i=0,1,...,n
  */
L1:
  fprintf(stdout, "Parameters for plane waves: a, n > ");
  scanf("%lf %d",&a, &n);
  if ( n < 1 || a <= 0) {
    fprintf(stderr, "wrong input parameters\n");
    exit(1);
  }
  fprintf (stdout,"   a, n = %f %d\n",a,n);
  /*
    Assign values of k_n
  */
  npw = 2*n+1;
  h  = (double *) malloc ( npw * npw * sizeof (double) );
  kn = (double *) malloc ( npw * sizeof (double) );
  e  = (double *) malloc ( npw * sizeof (double) );
  wrk= (double *) malloc (3*npw* sizeof (double) );

  kn[0] = 0.0;
  for ( i = 1; i < npw; i+=2 ) {
    kn[i  ] = (i+1)*pi/a;
    kn[i+1] =-(i+1)*pi/a;
  }
  for ( ij = 0; ij < npw*npw; ++ij ) {
    h[ij] = 0.0;
  }
  /*
    Assign values of the matrix elements of the hamiltonian 
    on the plane wave basis
  */
 
  ij = 0;
  for ( j = 0; j < npw; ++j ) {
    for ( i = 0; i < npw; ++i ) {
      /* NOTA BENE: the matrix h is a vector in the calling program,
         while dsyev expects a (pointer to a) fortran matrix.
         A fortran matrix is "simulated" in the following way:
         if h[ij] is the C array and  h(i,j) is a N*M Fortran matrix,
         h[ij] == h(i,j), where  ij = (j-1)*N + (i-1)  (i=1,N, j=1,M)
      */
      if ( i == j ) {
	h[ij++] = kn[i]*kn[i] - v0/a*b; 
      } else {
	h[ij++] = -v0/a * sin( (kn[j]-kn[i])*b/2.0 ) / (kn[j]-kn[i])*2.0;
      }
    }
  }
  /*
    Solution [expansion coefficients are stored into h(j,i)
    j=basis function index, i= eigenvalue index]
    (beware fortran-C reversed index convention!)
  */
  lwork = 3*npw;
  /* The leading dimension of array h is its first dimension, in
     this case npw, because columns are in consecutive order */
  dsyev_ ( V, U, &npw, h, &npw, e, wrk, &lwork, &info );
  if ( info != 0) {
    fprintf(stderr, "H-matrix diagonalization failed\n");
    exit(1);
  }
  printf("Lowest eigenvalues: %f %f %f\n",e[0],e[1],e[2]);
  /*
    Write to output file the lowest-energy state:
  */
  out = fopen("gs-wfc.out", "w");
  dx = 0.01;
  nr = (int) (a/2.0/dx+0.5);
  norm = 0.0;
  for ( i=-nr; i <= nr; i++ ) {
    x = dx*i;
    fr = 0.0; fi=0.0;
    /* the elements h[0]-h[npw-1] correspond to the first eigenvector */
    for ( j=0; j < npw; j++ ) {
      fr += h[j]*cos(kn[j]*x)/sqrt(a);
      fi += h[j]*sin(kn[j]*x)/sqrt(a);
    }
    prob = fr*fr+fi*fi;
    norm += prob*dx;
    fprintf(out,"%12.6f %10.6f %10.6f %10.6f\n",x, prob, fr, fi);
  }
  /*  verify normalization (if desired) */
  printf ("norm = %12.6f\n",norm);
  free(wrk); free(e); free(kn); free (h);
  goto L1;
}
