/* h2_hf_gauss */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define ABS(a)     (((a) < 0) ? -(a) : (a))
/* map fortran indices for matrices dimensioned ngaus into 1d arrays */
#define INDEX2(n,m) m*ngaus+n
#define INDEX4(i,j,n,m) m*pow(ngaus,3)+n*pow(ngaus,2)+j*ngaus+i

static double e2 = 2.0;
static double pi = 3.14159265358979;
static double tpi= 2.0*3.14159265358979;

main()
{
  /*
  
     Hartree-Fock ground-state solution for the H_2 molecule
     S=0 state - equivalent to Hartree approximation
     expansion on a gaussian basis set and self-consistency
     with diagonalization of the Fock matrix
     Ry atomic units: \hbar^2/2m = 1
     Requires lapack and blas
 */
  /* subroutines
  extern void diag (int , double *, double *, double *, double *); */

  /* variables */
  double zeta=1.0;
  unsigned int do_scf = 1;
  /*  do_scf=1 for H2 calculation, do_scf=0 for H2+ calculation */
  int ia,ib,ic,id, i,j,n, n_alpha, ngaus, iter, maxter = 100 ;
  int iab, iac, ibd, iabcd, iacbd ;
  double  aa, ab, dab, rab, d2, dmin, dmax, dd, deltae, eold, enew, enuc, e0 ;
  double *alpha, *d, *e, *c;
  double *h, *f, *s, *v, *q;
  char filename[80] = "h2.out";
  FILE *in, *out;

  fprintf(stdout, "Output to file %s\n",filename);
  if ( (out = fopen(filename, "w")) == NULL ) { 
     fprintf(stdout, "Error opening file\n");
     exit(1);
  }
  /*
        Input data
  */
  fprintf(stdout, "Parameters of the Gaussians from file >> ");
  scanf("%80s",filename); 
  if ( (in = fopen(filename, "r")) == NULL ) { 
     fprintf(stdout, "Error opening file\n");
     exit(1);
  }
  fprintf(stdout, "Reading from file %s\n",filename);
  n = fscanf(in,"%d", &n_alpha);
  if ( n_alpha < 1)  {
    fprintf(stderr, "n_alpha < 1, stopping\n");
    exit (0);
  } else {
    ngaus = 2*n_alpha;
    alpha =  (double *) malloc ( ngaus * sizeof (double) );
  }
  /* Read parameters of the gaussians */
  for ( ia = 0 ; ia < n_alpha ; ++ia ) {
    fscanf(in,"%lf", &alpha[ia]);
    fprintf (stdout,"Gaussian # %3d  coefficient >> %lf\n",ia, alpha[ia]);
  }
  fclose(in);
  /* Read distances */
  fprintf(stdout,"H2 distance (a.u.): Min, Max, step >> ");
  scanf("%lf %lf %lf", &dmin, &dmax, &dd);
  /*
     Starting solution (very lousy guess) for the first calculation
     Following calculations start from the solution at previous distance
  */
  c =  (double *) malloc ( ngaus * sizeof (double) );
  d =  (double *) malloc ( ngaus * sizeof (double) );
  e =  (double *) malloc ( ngaus * sizeof (double) );
  f =  (double *) malloc ( pow(ngaus,2) * sizeof (double) );
  v =  (double *) malloc ( pow(ngaus,2) * sizeof (double) );
  h =  (double *) malloc ( pow(ngaus,2) * sizeof (double) );
  s =  (double *) malloc ( pow(ngaus,2) * sizeof (double) );
  q =  (double *) malloc ( pow(ngaus,4) * sizeof (double) );
 
  for ( ia = 0 ; ia < ngaus ; ++ia ) { 
     c[ia] = 0.0;
  }
  c[0] = 1.0;
  /*
    loop over H-H distances (d2)
  */
  n = ( int ) ( (dmax-dmin)/dd + 0.5 );
  for ( j = 0 ; j <= n ; ++j) {        
     d2 = dmin + j*dd ;
     /*
       Basis set: b_i(r) = exp (-alpha_i (r-d_i)^2)
       (first n_alpha centered at -R/2, second n_alpha centered at +R/2)
     */
    for ( i = 0 ; i < n_alpha ; ++i) {        
        alpha[n_alpha+i] = alpha[i];
        d[i        ] = -d2/2.0;
        d[i+n_alpha] = +d2/2.0;
     }
     /*
           Assign values of overlap integrals S and of matrix elements 
           of the one-electron hamiltonian H on the gaussian basis:
     */
    for ( ia = 0 ; ia < ngaus ; ++ia ) {
       for ( ib = 0 ; ib < ngaus ; ++ib ) {
           aa = alpha[ia]+alpha[ib];
           ab = alpha[ia]*alpha[ib]/aa;
           dab= pow( d[ia]-d[ib], 2);
           /*  overlap */
           iab = INDEX2(ia,ib) ;
           s[iab] = pow( pi/aa, 1.5) * exp(-ab*dab);
           /* kinetic term of the hamiltonian */
           h[iab] = s[iab]*ab*(6.0-4.0*ab*dab);
           /* add coulombian term: beware the limit R_i-R_j=0 */
           for ( i=-1 ; i<=1; i=i+2 ) {
              rab= ABS ( ( alpha[ia]*d[ia]+alpha[ib]*d[ib] ) / aa - i*d2/2.0 );
              if ( rab < 1.0E-8 ) {
                 h[iab] = h[iab] - e2*zeta*tpi/aa;
              }
              else {
                 h[iab] = h[iab] - e2*zeta*s[iab] * erf( sqrt(aa)*rab ) / rab;
              }
           }
        }
     }
     /*
       Fill the Q matrix with the matrix elements of the e-e interaction:
       q(i,j,k,l) = \int b_i(r) b_j(r') 1/|r-r'| b_k(r) b_l(r') dr dr'
       beware the ordering of indices!
     */
    for ( ia = 0 ; ia < ngaus ; ++ia ) {
       for ( ib = 0 ; ib < ngaus ; ++ib ) {
          for ( ic = 0 ; ic < ngaus ; ++ic ) {
             for ( id = 0 ; id < ngaus ; ++id ) {
                 aa = alpha[ia]+alpha[ic];
                 dab= (alpha[ia]*d[ia]+alpha[ic]*d[ic])/aa;
                 ab = alpha[ib]+alpha[id];
                 rab= (alpha[ib]*d[ib]+alpha[id]*d[id])/ab;
                 rab= ABS ( dab-rab ) ;
                 /* iabcd is equivalent to (ia,ib,ic,id) in fortran */
                 iabcd = INDEX4(ia,ib,ic,id) ;
                 iac= INDEX2(ia,ic);
                 ibd= INDEX2(ib,id);
                 if ( rab < 1.0E-8 ) {
                    q[iabcd] = e2 * s[iac] * s[ibd] *
                         2.0 / sqrt(pi) * sqrt(aa*ab/(aa+ab)) ;
                 } else {
                    q[iabcd] = e2 * s[iac] * s[ibd] / rab *
                         erf ( sqrt(aa*ab/(aa+ab) )* rab ) ;
                 }
              }
          }
        }
     }
     /*
            Self-consistency iteration
     */
     enew = 0.0;
     deltae=1.0;
     for ( iter=1; (iter < maxter) && (deltae > 1.0e-8) ; iter++ ) {
        /*
                Fill the Fock matrix
        */
       for ( ia = 0 ; ia < ngaus ; ++ia ) {
          for ( ib = 0 ; ib < ngaus ; ++ib ) {
             iab= INDEX2(ia,ib) ;
             f[iab] = h[iab];
             if ( do_scf == 1 ) {
                for ( ic = 0 ; ic < ngaus ; ++ic ) {
                   for ( id = 0 ; id < ngaus ; ++id ) {
                      iacbd = INDEX4(ia,ic,ib,id) ;
                      /* the following is for Hartree only
                      f[iab] = f[iab] +  q[iacbd]*c[ic]*c[id];
                      */
                      /* the following is for Hartree-Fock */
                      iabcd = INDEX4(ia,ib,ic,id) ;
                      f[iab] = f[iab] + ( 2.0*q[iacbd] -
                                   q[iabcd] ) * c[ic]*c[id];
                   }
                }
             }
          }
       }
       /*
                Solution [expansion coefficients are stored into v(j,i)
                          j=basis function index, i= eigenvalue index]
       */
       diag ( ngaus, f, s,  e, v ) ;
        
       for ( i=0; i < ngaus; ++i) { 
          c[i] =  v[i];
       }
       eold = enew;
       enew = 0.0;
       for ( ia = 0; ia < ngaus; ++ia ) {
          for ( ib = 0; ib < ngaus; ++ib ) {
             iab= INDEX2(ia,ib) ;
             enew = enew + 2.0*h[iab]*c[ia]*c[ib];
             for ( ic = 0; ic < ngaus; ++ic ) {
                for ( id = 0; id < ngaus; ++id ) {
                      iacbd = INDEX4(ia,ic,ib,id) ;
                      /* the following is for Hartree only
                      enew = enew + q[iacbd]*c[ia]*c[ib]*c[ic]*c[id] ;
                      */
                      /* the following is for Hartree-Fock
                      */
                      iabcd = INDEX4(ia,ib,ic,id) ;
                      enew = enew + c[ia]*c[ib]*c[ic]*c[id] * 
                         ( 2.0*q[iacbd] - q[iabcd] ) ;
                }
             }
          }
       }
       /*
       fprintf (stdout, "Iteration # %d:  HF eigenvalue, energy = %12.6f %12.6f\n",
             iter, e[0], enew);
       */
       deltae = ABS ( enew - eold ) ;
       if ( deltae < 1.0e-8 ) {
          /* fprintf (stdout,"Convergence achieved, stopping\n"); */
          if ( j == 0 ) {
             fprintf (stdout,"     d      electron,    nuclear,    total (Ry) and binding (eV) energy\n");
          }
          if ( d2 > 1.0e-8 ) {
             enuc =  zeta*e2/d2;
          } else {
             enuc = 0.0;
          }
          if ( do_scf == 1 ) {
             e0   =-2.0;
          } else {
             e0   =-1.0; enew = e[0];
          }
          fprintf (stdout,"%8.6f %12.6f %12.6f %12.6f %12.6f\n",
              d2, enew, enuc, enew+enuc, (enew+enuc-e0)*13.6058);
          fprintf (out, "%8.6f %12.6f %12.6f %12.6f %12.6f\n",
              d2, enew, enuc, enew+enuc, (enew+enuc-e0)*13.6058);
       }
     }
  }
  fclose(out);
  free(q); free(s); free(h); free(v); free(f);
  free(e); free(d); free(c); free(alpha);
}

/* subroutine diag */
diag (int n, double *h, double *s, double *e, double *v)
{
  /*    Finds eigenvalues and eigenvectors of the generalized problem
	Hv=eSv, where H=hermitian matrix, S=overlap matrix */

  /* On input: n = dimension of the matrix to be diagonalized
               h = matrix to be diagonalized
               s = overlap matrix
     On output:
               e = eigenvalues
               v = eigenvectors
               s and h are unchanged */

  /* LOCAL variables */
  int lwork, i, j, k, nn, ij, info;
  /* lwork = dimension of workspace for lapack routine  */
  static double small = 1.0e-10;
  static char *V = "V";
  static char *U = "U";
  double *work, *aux, *uno;

  lwork=3*n;
  work = (double *) malloc( lwork * sizeof (double));

  /* Copy S into an auxiliary matrix (dsyev destroys the matrix) */
  aux = (double *) malloc( n * n * sizeof (double));
  for ( ij = 0; ij < n*n; ++ij ) {
    aux[ij] = s[ij];
  }

  /*  Diagonalize S  */

  dsyev_ ( V, U, &n, aux, &n, e, work, &lwork, &info ) ;

  if ( info !=0 ) {
    fprintf (stderr, "S-matrix diagonalization failed\n");
    exit (1);
  }

  /*    Keep only linearly independent combinations
	(within a given threshold)  */

  nn = 0;
  for ( i = 0; i < n; ++i ) {
    /*  i runs on all eigenstates
       nn runs on eigenstate with nonzero eigenvalue */
    if ( e[i] > small) {
      for ( j = 0; j < n; ++j ) {
	aux[j+nn*n] = aux[j+i*n] / sqrt(e[i]);
      }
      ++nn;
    }
  }

  if ( nn < n ) {
     fprintf (stdout, " # of linearly independent vectors = %d\n", nn);
   } 
  /*       Trasform H using the "aux" matrix
           V(i,j) = \sum_{k=1}^{n} H(i,k) aux(k,j),  i=1,n, j=1,nn
   */
  
  for ( i = 0; i < n; ++i ) {
    for ( j = 0; j < nn; ++j ) {
      v[i+j*n] = 0.0;
      for  ( k = 0; k < n; ++k ) {
	v[i+j*n] = v[i+j*n] + h[i+k*n] * aux[k+j*n];
      }
    }
  }
   /* h1(i,j) = \sum_{k=1}^{n} aux(k,i) v(k,j),  i=1,nn, j=1,nn
      H' = transpose(aux)*H*aux
   */
  uno = (double *) malloc( nn * nn * sizeof (double));
  for ( i = 0; i < nn; ++i ) {
    for ( j = 0; j < nn; ++j ) {
      uno[i+j*nn] = 0.0;
      for  ( k = 0; k < n; ++k ) {
	uno[i+j*nn] = uno[i+j*nn] + aux[k+i*n] * v[k+j*n];
      }
    }
  }
  
  /*  Diagonalize transformed H  */
  
  dsyev_ ("V", "U", &nn, uno, &nn, e, work, &lwork, &info );

  if ( info !=0 ) {
    fprintf (stderr, "H-matrix diagonalization failed\n");
    exit (1);
  }

  /*  Back-transform eigenvectors  */

  for ( i = 0; i < n; ++i ) {
    for ( j = 0; j < nn; ++j ) {
      v[i+j*n] = 0.0;
      for  ( k = 0; k < nn; ++k ) {
	v[i+j*n] = v[i+j*n] + aux[i+k*n] * uno[k+j*nn];
      }
    }
  }
  free(uno); free(aux); free(work);
  return;
}
