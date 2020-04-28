/* helium_gauss.c
   requires lapack
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define ABS(a)     (((a) < 0) ? -(a) : (a))
#define MIN(a,b)   (((a) < (b)) ? (a) : (b))
#define MAX(a,b)   (((a) > (b)) ? (a) : (b))

main()
{
  /*
    Hartree-Fock ground-state solution for the Helium atom
    S=0 state - equivalent to Hartree approximation
    expansion on a gaussian basis set and self-consistency
    with diagonalization of the Fock matrix
  */
  /* subroutine
  extern void diag (int , double *, double *, double *, double *); */
  static double e2 = 2.0;
  static double pi = 3.14159265358979;
  /* variables */
  int ia,ib,ic,id, i, n_alpha, iter, maxter=100;
  int iabcd, iacbd, iab;
  double zeta=2.0, aa, eold, enew;
  double *alpha, *e, *c;
  double *h, *f, *s, *v, *q;
  char filin[80];
  FILE *in, *out;

  /*   Input data */
  
  fprintf(stdout, " Atomic charge = %8.4f\n", zeta );
  fprintf(stdout, " Parameters of the Gaussians from file >> ");
  scanf("%80s", filin);
  if ( (in = fopen(filin, "r")) == NULL ) { 
     fprintf(stdout, "Error opening file\n");
     exit(1);
  }
  fscanf(in,"%d", &n_alpha);
  if ( n_alpha < 1)  {
    fprintf(stderr, "n_alpha < 1, stopping\n");
    exit (0);
  } else {
    alpha =  (double *) malloc ( n_alpha * sizeof (double) );
  }

  /* Read parameters of the gaussians */

  for ( ia = 0 ; ia < n_alpha ; ++ia ) {
    fscanf(in,"%lf", &alpha[ia]);
    fprintf (stdout," Gaussian # %3d  coefficient >> %lf\n",ia, alpha[ia]);
  }
  fclose(in);
  
  /*    Fill the Q matrix with the matrix elements of the e-e interaction:
        q(i,j,k,l) = \int b_i(r) b_j(r') 1/|r-r'| b_k(r) b_l(r') dr dr'
        where b_i(r) = exp (-alpha_i r^2) are the gaussians  */
  
  q =  (double *) malloc ( pow(n_alpha,4) * sizeof (double) );
  for ( ia = 0 ; ia < n_alpha ; ++ia ) {
    for ( ib = 0 ; ib < n_alpha ; ++ib ) {
      for ( ic = 0 ; ic < n_alpha ; ++ic ) {
        for ( id = 0 ; id < n_alpha ; ++id ) {
              /* iabcd is equivalent to (ia,ib,ic,id) in fortran */
              iabcd = id*pow(n_alpha,3)+ic*pow(n_alpha,2)+ib*n_alpha+ia;
              q[iabcd] = 4.0 * pow ( pi, 2.5 ) / 
                   ( alpha[ia]+alpha[ic])/(alpha[ib]+alpha[id] ) / 
                    sqrt( alpha[ia]+alpha[ib]+alpha[ic]+alpha[id] ) ;
           }
        }
     }
  }
  
  /*      Assign values of overlap integrals S and of matrix elements 
          of the one-electron hamiltonian H on the gaussian basis  */
  
  h =  (double *) malloc ( pow(n_alpha,2) * sizeof (double) );
  s =  (double *) malloc ( pow(n_alpha,2) * sizeof (double) );
  for ( ia = 0 ; ia < n_alpha ; ++ia ) {
    for ( ib = 0 ; ib < n_alpha ; ++ib ) {
      aa = alpha[ia] + alpha[ib];
      /* iab is equivalent to (ia,ib) in fortran */
      iab= ib*n_alpha+ia;
      s[iab] = pow( (pi/aa), 1.5 );
      h[iab] = s[iab]*6.0*alpha[ia]*alpha[ib]/aa - e2*zeta*2.*pi/aa;
    }
  }
  free (alpha);
 
  /*       Starting solution (very lousy guess)   */

  e =  (double *) malloc ( n_alpha * sizeof (double) );
  c =  (double *) malloc ( n_alpha * sizeof (double) );
  f =  (double *) malloc ( pow(n_alpha,2) * sizeof (double) );
  v =  (double *) malloc ( pow(n_alpha,2) * sizeof (double) );
  for ( ia = 1 ; ia < n_alpha ; ++ia ) {
     c[ia] = 0.0;
     }
  c[0] = 1.0;
  
  /*      Self-consistency iteration  */
  
  enew = 0.0;
  fprintf (stdout, "\n"); 
  for (iter = 1; iter < maxter ; ++iter){

     /*       Fill the Fock matrix  */
     
    for ( ia = 0 ; ia < n_alpha ; ++ia ) {
      for ( ib = 0 ; ib < n_alpha ; ++ib ) {
        iab= ib*n_alpha+ia;
        f[iab] = h[iab];
        for ( ic = 0 ; ic < n_alpha ; ++ic ) {
          for ( id = 0 ; id < n_alpha ; ++id ) {
            /* iacbd is equivalent to (ia,ic,ib,id) in fortran */
            iacbd = id*pow(n_alpha,3)+ib*pow(n_alpha,2)+ic*n_alpha+ia;
            /* the following is for Hartree only
            f[iab] = f[iab] + q[iacbd] * c[ic] * c[id] ;
            /* the following is for Hartree-Fock */
            /* iabcd = id*pow(n_alpha,3)+ic*pow(n_alpha,2)+ib*n_alpha+ia;
               f[iab] = f[iab] + ( 2.0*q[iacbd] - q[iabcd] ) * c[ic]*c[id]; */
          }
        }
      }
    }
    /*    Solution [expansion coefficients are stored into v(j,i)
                   j=basis function index, i= eigenvalue index]  */
  
    diag ( n_alpha, f, s,  e, v );
   
    for ( ia = 0 ; ia < n_alpha ; ++ia ) {
      c[ia] = v[ia];
    }
    eold = enew;
    enew = 0.0;
    for ( ia = 0 ; ia < n_alpha ; ++ia ) {
      for ( ib = 0 ; ib < n_alpha ; ++ib ) {
        iab= ib*n_alpha+ia;
        enew = enew + 2.0 * h[iab] * c[ia] * c[ib];
        for ( ic = 0 ; ic < n_alpha ; ++ic ) {
          for ( id = 0 ; id < n_alpha ; ++id ) {
            iabcd = id*pow(n_alpha,3)+ic*pow(n_alpha,2)+ib*n_alpha+ia;
            enew = enew +  q[iabcd] * c[ia] * c[ib] * c[ic] * c[id] ;
          }
        }
      }
    }
    fprintf (stdout, " Iteration # %2d: HF eigenvalue, energy: %f %f\n", 
       iter, e[0], enew );

    if ( ABS (enew-eold) < 1.0e-8){
        fprintf (stdout, "\n Convergence achieved, stopping\n");
        free(v); free(f); free(c); free(e); free(s); free(h); free(q);
        exit(0); 
     }
  }
  fprintf (stdout, "\n Convergence not reached, stopping\n");
  exit(2);
} /* end of main */

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
