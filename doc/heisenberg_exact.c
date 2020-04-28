
/* program heisenberg_exact */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#define MIN(a,b)   (((a) < (b)) ? (a) : (b))
#define INDEX(i,j) (j*nhil+i)

main()
{

/*
    Exact solution of isotropic 1-d Heisenberg model for s=1/2
    Periodic Boundary Conditions. Energy in units of J.
    Lanczos algorithm vs conventional Hamiltonian diagonalization
    (Hamiltonian stored as a matrix - limited to N<15 spins or so)
*/ 
  int Jsign = -1;  /* sign of J: -1 AFM, +1 FM  */
  unsigned int iseed = (unsigned int) time(NULL);
  /* seed for random number generation */
  int N, nup, nhil, nl, ncount, nonzero ;
  int ii, jj, lwork, info, i, j, ibip, one = 1 ;
  unsigned int k, *states;
  int *conf, *neigh ;
  double beta, uno = 1.0;
  double *H, *W, *work, *fat ;
  double *e, *d, *eaux, *daux, *v0, *v1, *v2 ;
  extern double dnrm2_ ( int *, double *, int * ) ;
  extern double ddot_ ( int *, double *, int *, double *, int * ) ;
  unsigned int search ( unsigned int *, int, unsigned int );
  !
  !
  fprintf (stdout,"Number of sites, number of up spins > ");
  scanf ("%d %d", &N,&nup) ;
  if ( N <= 0 ) { fprintf (stderr,"wrong N\n"); exit(1); }
  if ( N > 31 ) {
       fprintf (stderr,"N too big, not enough bits in an integer \n"); 
       exit(1);
  }
  /*
    fat(n) = n! - Should be integer but would overflow for N > 12
  */
  fat = (double *) malloc ( (N+1) * sizeof (double) );
  fat[0]=1 ;
  for ( i = 0 ; i < N ; i++ ) {
     fat[i+1]=fat[i]*(i+1) ;
  }
  /*
    nhil = N!/nup!/ndw! is the dimension of the Hilbert space
    Not a smart implementation: nhil is integer by construction
  */
  nhil = ( int ) ( fat[N] /(fat[nup]*fat[N-nup]) + 0.5 ) ;
  free ( fat ) ;
  
  conf = (int *) malloc ( N * sizeof (int) );
  neigh= (int *) malloc ( N * sizeof (int) );
  states=(int *) malloc ( nhil * sizeof (int) );
  H    = (double *) malloc ( nhil * nhil * sizeof (double) );

  for ( j = 0; j < N; j++) {
     /*  neigh(j) is the index of the neighbor of the j-th spin
                  (only the one at the right to prevent double counting) */
     neigh[j] = j+1 ;
  }
  /*  periodic boundary conditions */
  neigh[N-1] = 0 ;

  ncount=0;
  /* the do loop runs on all possible 2^N states
     the bit representation is used to represent up and down spins
  */
  for ( k = 0; k < pow(2,N); k++ ) {
     /*  count all up spin  */
      ibip = 0 ;
      for ( j = 0; j < N; j++ ) {
         /*  btest(k,j) = .true. if bit j in k is 1 (i.e. spin j is up) */
         if ( btest(k,j) ) { ibip++ ; }
      }
      if ( ibip == nup ) {
         /* k labels a basis set state with Nup spin up */
         /* fprintf(stdout,"state(%d)=%d\n",ncount,k); */
         states[ncount++] = k ;
      }
      if ( ncount > nhil ) {
         fprintf ( stderr, "Dimension of Hilbert space: %d, is insufficient\n",
                   nhil ) ;
         exit(1);
      }
   }
   if ( ncount != nhil ) {
      fprintf ( stderr, "Dimension of Hilbert space mismatch: %d %d\n", 
                nhil, ncount ) ;
      exit(1);
   } else {
      fprintf ( stdout, "Dimension of Hilbert space: %d\n", nhil) ;
   }

   for ( i = 0; i < nhil*nhil; i++ ) { H[i] = 0.0; }

   /* Fill the Hamiltonian: (J/2)\sum_{ij} S_-(i)S_+(j) + h.c. term */
   nonzero = 0 ;
   for ( ii = 0; ii < nhil; ii++ ) {
      /*  do on all states (ii) of the Hilbert space */
      for ( j = 0;  j < N; j++ ) {
         /* conf(j) = -1 if j-th bit of state ii is 0 (i.e. spin down)
            conf(j) = +1 if j-th bit of state ii is 1 (i.e. spin up)  */
         conf[j] = -1;
         if ( btest(states[ii],j) ) { conf[j] = 1 ; }
         /* fprintf(stdout,"bit %d of state %d: %d\n",j,ii,conf[j]); */
      } 
      /*  conf(j) = sequence of + and - spins for state ii */
      for ( j=0; j < N; j++ ) {
         i=neigh[j];
         /* i is the index of the neighbor of spin j-th */
         if ( conf[j] == -1 && conf[i] == 1 ) {
            /* state |...sigma_i...sigma_j...> with sigma_j=-1, sigma_i=+1
               nonzero contribution from term (1/2)\sum_{ij} S_-(i)S_+(j) */
            k = states[ii] + pow(2,j)-pow(2,i) ;
            /*  k is the label of state <...sigma_i...sigma_j...|
                                  with sigma_j=+1, sigma_i=-1  */
            jj = search ( states, nhil, k ) ;
            /* jj is the index of state k: states(jj) = k */
            H[ INDEX(ii,jj) ] = H[ INDEX(ii,jj) ] - 0.5*Jsign ;
            nonzero++ ;
         } else if ( conf[j] == 1 && conf[i] == -1 ) {
            /*  as above for (1/2)\sum_{ij} S_-(i)S_+(j)  */
            k = states[ii] - pow(2,j) + pow(2,i) ;
            jj = search ( states, nhil, k ) ;
            H[ INDEX(ii,jj) ] = H[ INDEX(ii,jj) ] - 0.5*Jsign ;
            nonzero++;
          }
       }
    }
    /*  fill the Hamiltonian: J\sum_{ij} S_z(i)S_z(j) term
        (contributes only to the diagonal of H)             */
    for ( ii = 0; ii < nhil; ii++) {
       for ( j = 0; j < N; j++ ) {
          conf[j] = -1;
          if ( btest(states[ii],j) ) conf[j]=1 ;
       }
       for ( j = 0 ; j < N; j++ ) {
          i=neigh[j];
          H[ INDEX(ii,ii) ] = H[ INDEX(ii,ii) ] - Jsign*0.25*conf[j]*conf[i] ;
       }
       nonzero++;
    }

    fprintf(stdout,"Number of nonzero elements: %d (%6.1f %% of total)\n",
          nonzero, (double) (nonzero)/nhil/nhil*100 ) ;
    
    fprintf(stdout,"Number of Lanczos steps > ");
    scanf ("%d", &nl) ;
    nl = MIN ( nl, nhil ) ;
    /* d, e vectors containing the Hamiltonian in tridiagonal form */
    d = (double *) malloc ( (nl+1) * sizeof (double) );
    e = (double *) malloc ( (nl+1) * sizeof (double) );
    /* vectors used in the Lanczos algorithm */
    v0 = (double *) malloc ( nhil * sizeof (double) );
    v1 = (double *) malloc ( nhil * sizeof (double) );
    v2 = (double *) malloc ( nhil * sizeof (double) );
    /* fill starting vector with random numbers */
    srand (iseed);
    for ( i=0; i < nhil; i++ ) { v1[i] = rand() ; }
    e[0] = dnrm2_ ( &nhil, v1, &one );
    for ( i=0; i < nhil; i++ ) { v1[i] = v1[i] / e[0]; v0[i] = 0.0;  }
    e[0] = 0.0; d[0] = 0.0 ;
    /* Lanczos procedure starts here  */
    for ( j = 1; j <= nl ; j++ ) {
       /*  v_{j-1} == v0
           v_j     == v1
           v_{j+1} = H*v_j - beta*v_{j-1}  (overwritten to v0)
       */ 
       beta = -e[j-1];
       dgemv_ ( "N", &nhil, &nhil, &uno, H, &nhil, v1, &one, &beta, v0, &one ) ;
       /*  alpha = <v_{j+1}|v_j>  */
       d[j] = ddot_ ( &nhil, v0, &one, v1, &one ) ;
       /*  v_{j+1} = v_{j+1} - alpha*v_j   (is orthogonal to v_j) */ 
       for ( i=0; i < nhil; i++ ) { v2[i] = v0[i] - d[j]*v1[i]; }
       /* beta = |v_{j+1}| */
       e[j] = dnrm2_ ( &nhil, v2, &one )  ;
       /*  v_{j+1} = v_j, v_j = v_{j+1}/beta  */
       for ( i=0; i < nhil; i++ ) { v0[i] = v1[i]; v1[i] = v2[i] / e[j]; } 
    }
    free (v2); free(v1); free(v0);

    /*  Hamiltonian is now in tridiagonal form: solve */

    daux = (double *) malloc ( nl * sizeof (double) );
    eaux = (double *) malloc ( nl * sizeof (double) );
    
    for ( j=1; j <= nl; j++ ) {
       /*  arrays containing diagonal and subdiagonal are destroyed by dsterf */
       for ( i=1; i <= nl; i++ ) { daux[i-1] = d[i]; eaux[i-1]=e[i] ; }
       /*  note that d(0) and e(0) are zero and are not copied  */
       /*  dsterf calculates eigenvalues only of a tridiagonal matrix  */
       dsterf_ (&j, daux, eaux, &info ) ;
       if (info != 0) {
          fprintf(stdout,"Error in dsterf, info=%d\n", info);
       }
       fprintf(stdout,"nl=%d,  E=%lf,   E/N=%lf\n", j,daux[0],daux[0]/N );
    }
    free (eaux) ; free (daux); free (e); free(d);

    /*  check: conventional diagonalization */

    lwork = 3*nhil;
    W = (double *) malloc ( nhil * sizeof (double) );
    work = (double *) malloc ( lwork * sizeof (double) );

    dsyev_ ("N","U", &nhil, H, &nhil, W, work, &lwork, &info ) ;
    free (work) ;
    if (info != 0) {
          fprintf(stdout,"Error in dsyev, info=%d\n", info);
    }
    fprintf(stdout,"exact:  E=%lf,   E/N=%lf\n",W[0],W[0]/N );
    
    free (W); free (states); free (H); free (neigh); free (conf);
  }

unsigned int search ( unsigned int* p, int nhil, unsigned int k )  {
    /*
      On input:
           p(1:nhil) = labels for all states (ordered: p(i+1) > p(i) )
           k         = a state label
      On output:
           index i such that p(i) = k
    unsigned int *p;
    int nhil;
    unsigned int k;
    */
    unsigned int imin, imax, lim, l, ii;

    lim = log( (double ) (nhil))/log(2.0) + 1 ;
    /* 2^lim > Nhil: lim = max number of steps needed to locate p(i)=k  */

    imin = 0;
    imax = nhil-1;
    for ( l = 0; l < lim; l++ ) {
       ii=(imin+imax)/2 ;
       if (p[ii] == k) {
          imin=ii ;
          return ii ;
       } else {
          if ( p[ii] > k ) {
             imax=ii-1;
          } else {
             imin=ii+1;
          }
       }
    }

    fprintf(stderr, "Something wrong: search not converged\n");
    exit(1);
}

int btest(unsigned int k , int j)
{
    unsigned int uno = 1;
    return (k >> j) & uno ;
}

