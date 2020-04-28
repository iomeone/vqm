/* program cohenbergstresser */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define ABS(a)     (((a) < 0) ? -(a) : (a))
#define MIN(a,b)   (((a) < (b)) ? (a) : (b))
#define MAX(a,b)   (((a) > (b)) ? (a) : (b))

static double e2 = 2.0;
static double pi = 3.14159265358979;
static double tpi= 2.0*3.14159265358979;

main()
{
  /*
        Band structure of semiconductors in the zincblende and
        diamond structure - Cohen and Bergstresser, PRB 141, 789 (1966)
        Expansion on a plane-wave basis set and diagonalization
        Units: hbar^2/1m = 1
        Requires lapack dsyev
  
        a  = lattice parameter (a.u.)
        V*N= pseudopotential form factors, in Ry, as in the paper
         s = symmetric term, a = antisymmetric term
            The following numbers are good for Silicon
  */
  double a = 10.26, vs3 = -0.21, vs8 = 0.04, vs11 = 0.08, 
                    va3 = -0.00, va4 = 0.00, va11 = 0.00 ;
  /*    the two atoms are at -tau and +tau (in units of a) */
  double tau1 = 0.125, tau2 = 0.125, tau3 = 0.125 ;
                      
  int n, npw, nk;
  double ecut, kg2, g2, vag, vsg, k[3], gx,gy,gz, h1[3], h2[3], h3[3];
  double *kx, *ky, *kz, *kgx, *kgy, *kgz, *e, *work;
  double *h;
  int i,j, ij, m, nmax, n1, n2, n3, lwork, info;
  FILE *in, *out;
  /*
       Plane waves basis set: G_n=n*2pi/L, \hbar^2/2m*G^2 < Ecut
  */
  fprintf (stdout, "Cutoff for plane waves: ecut (Ry) > ");
  scanf ("%lf", &ecut) ;
  if ( ecut <= 0.0 ) { fprintf (stderr," wrong cutoff\n"); exit(1); }
  /*
       Number and list of k-vectors
  */
  fprintf (stdout,"Number of k-vectors > ");
  scanf ("%d", &nk) ;
  if ( nk <= 0 ) { fprintf (stderr," wrong input parameter\n"); exit(1); }
  kx = (double *) malloc ( nk * sizeof (double) );
  ky = (double *) malloc ( nk * sizeof (double) );
  kz = (double *) malloc ( nk * sizeof (double) );
  fprintf (stdout,"k (in 2pi/a units) > \n");

  for ( i = 0 ; i < nk ; i++ ) {
     scanf ("%lf %lf %lf", &kx[i], &ky[i], &kz[i] ) ;
  }
  /*
       Basis vectors for reciprocal lattice in 2pi/a units
  */
  h1[0] = 1.0;   h1[1] = 1.0; h1[2] =-1.0;
  h2[0] = 1.0;   h2[1] =-1.0; h2[2] = 1.0;
  h3[0] =-1.0;   h3[1] = 1.0; h3[2] = 1.0;
  /*
       Loop on k-vectors
  */ 
  for ( n = 0 ; n < nk ; n++ ) {
     /*
         Count plane waves such that (\hbar^2/2m)(k+G)^2 < Ecut
         nmax is an estimate of max index useful in the generation
         of PWs as G(n1,n2,n3) =  n1*h1 + n2*h2 + n3*h3
     */
     nmax = (int) ( sqrt (ecut) / (tpi/a * sqrt(3.0) ) + 0.5 ) + 1 ;
     npw = 0;
     for ( n1 = -nmax; n1 <= nmax; n1++ ) {
        for ( n2 = -nmax; n2 <= nmax; n2++ ) {
           for ( n3 = -nmax; n3 <= nmax; n3++ ) {
              k[0] = kx[n] ; k[1] = ky[n] ; k[2] = kz[n] ; 
              kg2 = 0.0 ;
              for ( i = 0; i < 3 ; i++ ) { 
                 kg2 = kg2 + pow( k[i] + n1*h1[i] + n2*h2[i] + n3*h3[i], 2);
              }
              /*  kg2 is now |k+G|^2 in 2pi/a units  */
              kg2 = pow(tpi/a,2)*kg2 ; 
              if ( kg2 <= ecut ) { ++npw ; }
           }
        }
     }
     if ( npw < 1 ) {
        fprintf (stdout, "Incorrect number of plane waves (%d)!\n",npw);
        exit(1);
     } else {
        fprintf (stdout, "Number of plane waves=%d\n",npw);
     }
     h = (double *) malloc ( npw*npw*sizeof (double) );
     kgx = (double *) malloc ( npw * sizeof (double) );
     kgy = (double *) malloc ( npw * sizeof (double) );
     kgz = (double *) malloc ( npw * sizeof (double) );
     e   = (double *) malloc ( npw * sizeof (double) );
     lwork = 3*npw;
     work= (double *) malloc ( lwork*sizeof (double) );
     /*
         now generate PWs (beware: in 2pi/a units)
     */
     j = 0;
     for ( n1 = -nmax; n1 <= nmax; n1++ ) {
        for ( n2 = -nmax; n2 <= nmax; n2++ ) {
           for ( n3 = -nmax; n3 <= nmax; n3++ ) {
              k[0] = kx[n] ; k[1] = ky[n] ; k[2] = kz[n] ;
              for ( i = 0; i < 3 ; i++ ) { 
                 k[i] = k[i] + n1*h1[i] + n2*h2[i] + n3*h3[i] ;
              }
              kg2 = pow(tpi/a,2) * ( k[0]*k[0] + k[1]*k[1] + k[2]*k[2] );
              if ( kg2 <= ecut ) {
                 kgx[j] = k[0] ;  kgy[j] = k[1] ;  kgz[j] = k[2] ; j++ ;
              } 
           }
        }
     }
     if ( j != npw ) { 
        fprintf (stderr,"Some PWs are missing (%d < %d)\n",j+1,npw);
        exit(1);
     }
     /* cleanup */
     ij = 0 ;
     for ( i = 0 ; i < npw ; i++ ) {
        for ( j = 0 ; j < npw ; j++ ) { h[ij++]= 0.0 ; }
     }
     /*
            Assign values of the matrix elements of the hamiltonian 
            on the plane wave basis
     */
     ij = 0 ;
     for ( i = 0 ; i < npw ; i++ ) {
        for ( j = 0 ; j < npw ; j++ ) {
           gx = kgx[i] - kgx[j] ;
           gy = kgy[i] - kgy[j] ;
           gz = kgz[i] - kgz[j] ;
           g2 = gx*gx + gy*gy + gz*gz ;
           if ( ABS (g2-3.0) < 1.0e-6 ) {
              vsg = vs3;
              vag = va3;
           } else if ( ABS (g2-4.0) < 1.0e-6 ) {
              vsg = 0.0;
              vag = va4;
           } else if ( ABS (g2-8.0) < 1.0e-6 ) {
              vsg = vs8;
              vag = 0.0;
           } else if ( abs (g2-11.0) < 1.0e-6 ) {
              vsg = vs11;
              vag = va11;
           } else {
              vsg =0.0;
              vag =0.0;
           }
           if ( i == j ) {
              h[ij++] = vsg +  pow ( tpi/a, 2) *
                       ( pow (kgx[i],2) + pow (kgy[i],2) + pow (kgz[i],2) ) ;
           } else {
              h[ij++] = vsg * cos(tpi*(gx*tau1+gy*tau2+gz*tau3)) ;
                     /* vag * sin(tpi*(gx*tau1+gy*tau2+gz*tau3) ) */
           }
           /* print  '(2i4,f12.6)', i,j, h(i,j) */
        }
     }
     /*
            Solution [expansion coefficients are stored into h(j,i)
                      j=basis function index, i= eigenvalue index]
     */
     dsyev_ ( "V", "U", &npw, h, &npw, e, work, &lwork, &info );
     if ( info != 0 ) { 
        fprintf (stderr,"H-matrix diagonalization failed\n") ;
        exit(1);
     }
     for ( i = 0 ; i < npw ; i++ ) {
         e[i] = e[i]*13.6058;
     }
     fprintf(stdout, "k =  %7.4lf %7.4lf %7.4lf\n",kx[n],ky[n],kz[n]);
     fprintf(stdout, "%12.6lf %12.6lf %12.6lf %12.6lf\n",e[0],e[1],e[2],e[3]);
     fprintf(stdout, "%12.6lf %12.6lf %12.6lf %12.6lf\n",e[4],e[5],e[6],e[7]);
     /*
           Write to output file the band dispersion e(k)
     */
     free(h); free (work) ; free(e) ; free (kgz) ; free (kgy) ; free (kgx);
  }
  free (kz); free (ky); free (kx) ;
}
