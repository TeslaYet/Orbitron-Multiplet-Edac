#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "w3j.h"
#include "globals.h"
#include "ndimarraycalloc.h"

#define PI 3.14159265358979323844
#define SQRTHALF 0.707106781186547524

#define EPSILON 1.e-10
#define EPSPES 1.e-10

void printstatevector( int nstates, double *vec, struct Fock *state ) ;

int printspec( int dim, double *lambda, int *pgstdeg, double *pgstenergy ) ;


void thetaphimesh( double thmin, double thmax, double thdelta, 
                   double phmin, double phmax, double phdelta, 
                int *pntheta, double **ptheta, int **pnphi, double ***pphi ) ;
void printrpesmatele( FILE *fp, double omega, int npesorb, int nnm1fst,
		      int gstdeg, double *gstweight, double *nm1fenergy,
		      double gstenergy, complex ****rpesmatele ) ;
void printrpesxmat( FILE *fp, double omega, int npesorb, int nnm1fst,
		    int gstdeg, double *gstweight, complex ****rpesmatele ) ;
int main() {
  const int ncvsh = 2 ; /* BAD PROGRAMMING */
  int nshells, nconfs, nelectrons ;
  int i, j, k, npesorb, isorb, iso, jso, m, ksum = 0 ;
  int *lsh,  *sorb1sh, **occ;
  int ist, jst, fist, fjst, nstates, nnm1fst, nfstates ;
  int nhamele, degeneracy, info, gstdeg, nlevels, fst0deg ;
  int inishell, q, iq, nlistele, ngstbasis, nmstates ;
  int igstdeg, ntheta, *nphi, ngam, nomega ;
  double *ksish, gstenergy, **gstvec, *fstenergy, fst0energy, oldlambda ;
  double e2p, e3d, *gstweight ;
  double **cfdmat, hmag, thetamag, sum, dum, buf[3], *nt2g, *neg, *ndxy ;
  double *ham, *lambda, **nm1fstvec, *nm1fenergy, **mstvec, *menergy, *fstvec ;
  double ****pesmatele, ***xasmatele, **fbasvmst, **fstvmst ;
  double **radipmatele, omega, ommin, ommax, deltaom, gamma ; 
  double *egam, *gam, gamma0, *gamst ;
  double thmax, thmin, thdelta, phmax, phmin, phdelta, *theta, **phi ; 
  double sumr, suml;
  double eps[3][3]={{.5,SQRTHALF,.5},{-SQRTHALF,0.,SQRTHALF},{.5,-SQRTHALF,.5}};
  double complex csum, *cmdum, csumxaq[3]; 
  complex ****rpesmatele, zdum ;
  struct Fock *state, *gstbasis, *fstate, *nm1fstbas, *mstbas ;
  struct Spamaline *hamsparse, *hamsparsep, *pdipsmline0[3], *pvaismline0 ;
  struct Spamaline *pfvmsml0, *psml, *pt2gsmline0, *pegsmline0, *pdxysmline0 ;
  struct O1plistitem *po1plistitem ;	
  struct O2plistitem *po2plistitem ;	
  struct O1p dipop[3], t2gop, egop, dxyop ;
  struct O2p vaiop ;
  FILE *fp, *fpx, *fpy, *fpc, *fpa, *fpz ;


  w3jtabmake(); 
/*
  printf("Enter thmin, thmax, thdelta, phmin, phmax, phdelta\n") ;
  scanf("%lf%lf%lf%lf%lf%lf", 
         &thmin, &thmax, &thdelta, &phmin, &phmax, &phdelta) ; 
*/
/* make emission-angle mesh theta[i],phi[i][j], i<ntheta, j<nphi[i]
  thetaphimesh( thmin, thmax, thdelta, phmin, phmax, phdelta, 
                &ntheta, &theta, &nphi, &phi ) ;
  for ( i = 0 ; i < ntheta ; i++ ) {
    printf("%6.2lf, [%2d]", theta[i], nphi[i] ) ;
    for ( j = 0 ; j < nphi[i] ; j++ ) printf(" %6.2lf", phi[i][j] ) ;
    printf("\n");
  }
*/
  cfdmat = (double **) malloc( 5 * sizeof(double *) );
  for (i=0;i<5;i++) cfdmat[i] = (double *) malloc( 5 * sizeof(double) );
  printf("Enter e2p, e3d:") ; scanf("%lf%lf", &e2p, &e3d ) ;
  printf("Enter CF on d-shell as real 5x5 matrix:\n") ;
  for ( i = 0 ; i < 5 ; i++ )
    for ( j = 0 ; j < 5 ; j++ )
      scanf("%lf", &cfdmat[i][j] ) ;
  printf("Enter (mag.field) h, theta(in deg): ") ;
  scanf("%lf%lf", &hmag, &thetamag ) ;
  thetamag *= PI/180. ;
  
  printf("Enter ommin, ommax, deltaom, gamma0, ngam, egam[0], gam[0] .., gam[ngam-1] ") ;
  scanf("%lf%lf%lf%lf%d", &ommin, &ommax, &deltaom, &gamma0, &ngam ) ;
  nomega = (int) ( (ommax-ommin)/deltaom + 1.00001 ) ;  
  egam=(double *)malloc(ngam*sizeof(double));
  gam =(double *)malloc(ngam*sizeof(double));
  for ( i = 0 ; i < ngam ; i++ )  scanf("%lf%lf", &egam[i], &gam[i] ) ;
  readshells( &nshells, &lsh, &sorb1sh, &ksish) ;
  radipmatele = calloc2double( nshells, nshells) ;
  printf("Enter radipmatele's shells: 0->1, 1->2, 1->3: ") ;
  if ( nshells > 3 ) scanf("%lf%lf%lf",
      &radipmatele[0][1],&radipmatele[1][2],&radipmatele[1][3] );
  else { printf("radip trouble\n"); exit(1);}

  readconfs( nshells, lsh, &nconfs, &nelectrons, &occ, &nstates ) ;

  makestates( nshells, lsh, sorb1sh, nconfs, nelectrons, occ, nstates,
	      &state) ;
  nhamele = calcham( &hamsparse, nstates, state, nconfs, nshells, lsh, 
		     sorb1sh, occ, ksish, e2p, e3d, cfdmat, hmag, thetamag ) ;

  for ( i = 0 ; i < nconfs ; i++ )
    free( occ[i] ) ;
  free( occ ) ;

  printf("nhamele = %d\n", nhamele ) ;
  
  ham = ( double * ) calloc( nstates * nstates, sizeof( double ) ) ;
  lambda = ( double * ) malloc( nstates * sizeof( double ) ) ;
 
  for ( i = 0 ; i < nstates * nstates ; i++ )
    ham[i] = 0. ;

  hamsparsep = hamsparse ; 
  for ( ist = 0 ; ist < nstates ; ist++ ) {
    for ( k = 0 ; k < hamsparsep -> n ; k++ ) 
      ham[ ist + nstates * hamsparsep -> j[k] ] = hamsparsep -> v[k] ;
    hamsparsep++ ;
  }

  info =  solve_dsyev( nstates, ham, lambda ) ; 
  if ( info ) { printf(" DSYEV : INFO = %d\n", info ) ; exit(1) ; } 


   
  spamadelete( hamsparse, nstates ) ; 
  
  nlevels = printspec( nstates, lambda, &gstdeg, &gstenergy ) ; 

  printf("nlevels = %d\n", nlevels ) ;

  /* store degenerate ground states in gstvec[gstdeg][nstates] */

  gstweight = ( double * ) calloc( gstdeg, sizeof( double ) ) ;
  for ( k = 0 ; k < gstdeg ; k++ ) gstweight[k] = 1. / ( (double) gstdeg ) ;
  gstvec = ( double ** ) malloc( gstdeg * sizeof( double * ) ) ;
  for ( k = 0 ; k < gstdeg ; k++ ) {
    gstvec[k] = ( double * ) malloc( nstates * sizeof( double ) ) ; 
    for ( ist = 0 ; ist < nstates ; ist++ )
      gstvec[k][ist] = ham[ ist + nstates * k ] ;
  }
  printf("gstenergy = %lf.  %d degen. gstvec's:\n", gstenergy, gstdeg ) ;
  for ( k = 0 ; k < gstdeg ; k++ ) 
    printstatevector( nstates, gstvec[k], state ) ;

  ngstbasis = nstates ;
  gstbasis = state ;

  free( lambda ) ;
  free( ham ) ;
  

  readconfs( nshells, lsh, &nconfs, &nelectrons, &occ, &nstates ) ;

  makestates( nshells, lsh, sorb1sh, nconfs, nelectrons, occ, nstates,
	      &state) ;
  nhamele = calcham( &hamsparse, nstates, state, nconfs, nshells, lsh, 
		     sorb1sh, occ, ksish, e2p, e3d, cfdmat, hmag, thetamag ) ;
  printf("nhamele = %d\n", nhamele ) ;

  for ( i = 0 ; i < nconfs ; i++ )
    free( occ[i] ) ;
  free( occ ) ;

  printf("after free occ\n");


  ham = ( double * ) calloc( nstates * nstates, sizeof( double ) ) ;
  nm1fenergy = ( double * ) malloc( nstates * sizeof( double ) ) ;

  for ( i = 0 ; i < nstates * nstates ; i++ )
    ham[i] = 0. ;

  hamsparsep = hamsparse ; 
  for ( ist = 0 ; ist < nstates ; ist++ ) {
    for ( k = 0 ; k < hamsparsep -> n ; k++ ) 
      ham[ ist + nstates * hamsparsep -> j[k] ] = hamsparsep -> v[k] ;
    hamsparsep++ ;
  }
  spamadelete( hamsparse, nstates ) ; 

  info =  solve_dsyev( nstates, ham, nm1fenergy ) ; 
  if ( info ) { printf(" DSYEV : INFO = %d\n", info ) ; exit(1) ; } 

  nnm1fst = nstates ;
  nm1fstbas = state ;  /* keep this N-1 final state basis */

 /* store nm1-final state vectors */
  nm1fstvec = ( double ** ) malloc( nnm1fst * sizeof( double * ) ) ;
  for ( jst = 0 ; jst < nnm1fst ; jst++ ) 
    nm1fstvec[jst] = ham + nnm1fst * jst ;        

/* long way ... safer ? ...
  nm1fstvec = calloc2double(nnm1fst,nnm1fst) ;
  for ( jst = 0 ; jst < nnm1fst ; jst++ ) 
    for ( ist = 0 ; ist < nnm1fst ; ist++ ) {
      nm1fstvec[jst][ist] = ham[ist + nnm1fst * jst] ;        
    }
  free( ham ) ;
*/


  nlevels = printspec( nstates, nm1fenergy, &fst0deg, &fst0energy ) ; 

  printf("nlevels = %d\n", nlevels ) ;
  for ( i = 0 ; i < nstates ; i++ ) { 
    printf(" %5.3lf", nm1fenergy[i]) ; 
    if ( i%20==19 ) printf("\n"); 
  }
  printf("\n"); 
  
/* from |N-1> final states to |N-1>|PE> final states */
  npesorb = sorb1sh[nshells] - sorb1sh[ncvsh] ;
  nfstates = nnm1fst * npesorb ; 
  fstate = ( struct Fock * ) calloc( nfstates, sizeof( struct Fock ) ) ;
  fstenergy = ( double * ) calloc( nfstates, sizeof( double ) ) ;
  fstvec = ( double * ) malloc( nfstates * sizeof( double ) ) ;
  for ( iso = 0 ; iso < npesorb ; iso++ ) {
    isorb = sorb1sh[ncvsh] + iso  ;
    for ( ist = 0 ; ist < nnm1fst ; ist++ ) {
      fist = iso*nnm1fst + ist ;
      fstate[fist] = state[ist] ;
      fockcre( fstate + fist, isorb ) ;
      printf("%3d:%18llu)=", fist, fstate[fist].n) ; 
      fockdisplay( fstate + fist ) ; printf("\n") ;
/* Note: Energy of PE=   ep[fist] = gstenergy + omega - nm1fenergy[ist] ; */
      fstenergy[fist] = nm1fenergy[ist] ; 
    }
  }


/* construct dipole operator  P^(1)_q = r C^(1)_q , q=-1,0,1 [Cowan (14.24)] */
  inishell = 1 ;
  for ( iq = 0 ; iq < 3 ; iq++ ) {
    q = iq - 1 ;
    po1plistitem = 0 ;
    nlistele = dipole2list( nshells, lsh, sorb1sh, inishell, q, radipmatele,
			       &po1plistitem) ;
    printf("Dipole Op. q = %d. nlistele = %d.\n", q, nlistele ) ;
    dipop[iq] = o1pmake( po1plistitem ) ;
    o1plistdelete( po1plistitem ) ;
    nhamele = o1ptospama( ngstbasis, gstbasis, nfstates, fstate, 
			  dipop[iq], &pdipsmline0[iq] ) ;
    printf("numtotele = %d\n", nhamele ) ;
    for ( i = 0 ;  i < nfstates ; i++ ) {
      printf("<%3d|r_{%2d}| . > . :", i, q ) ;   
      spamalinewrite( pdipsmline0[iq][i] ) ;
    }
  }

  pesmatele = calloc4double( 3, npesorb, gstdeg, nnm1fst ) ;
  for( jso = 0 ; jso < npesorb ; jso++ ) {
    for ( jst = 0 ; jst < nnm1fst ; jst++ ) {
      fjst = jso*nnm1fst + jst ;
/* generally: eigenvector  |j)  = Sum_i |bas_i> C_ij, where C_ij=ham[i+n*j]
Here: direct product: |j)|jso) = Sum_{i,iso} |bas_i>|iso> C_{i,iso;j,jso}
                               = Sum_i |bas_i>|jso> C_ij
*/
      for ( fist = 0 ; fist < nfstates ; fist++ ) fstvec[fist] = 0. ;
      for ( ist = 0 ; ist < nnm1fst ; ist++ ) 
        fstvec[jso*nnm1fst + ist] = nm1fstvec[jst][ist] ; /* = ham[ist + nnm1fst * jst] */
      for ( iq = 0 ; iq < 3 ; iq++ ) {
        for ( igstdeg = 0 ; igstdeg < gstdeg ; igstdeg++ ) {
          pesmatele[iq][jso][igstdeg][jst] =
          spamarealmatele( ngstbasis, gstvec[igstdeg], nfstates, fstvec, 
  			 pdipsmline0[iq] ) ;
        }
      }
    }
  }
/* free pdipsmline0 */
  for ( iq = 0 ; iq < 3 ; iq++ ) 
    spamadelete( pdipsmline0[iq], nfstates ) ;

/*  print out total PES spectrum Sum_q Sum_lms */
  fp = fopen("pes.dat","w") ;
  for ( ist = 0 ; ist <  nnm1fst ; ist++ ) {
    sum = 0. ;
    for ( iq = 0 ; iq < 3 ; iq++ ) 
      for ( iso = 0 ; iso < npesorb ; iso++ ) 
        for ( igstdeg = 0 ; igstdeg < gstdeg ; igstdeg++ ) {
          dum = pesmatele[iq][iso][igstdeg][ist] ;
          sum += dum*dum ;
	}
    if ( sum > EPSPES) {
      fprintf(fp,"%11.6lf %15.8e\n",fstenergy[ist]-gstenergy, sum/gstdeg ) ;
     /* picks up fist = ist (iso=0) */
    }
  }   
  fclose(fp) ;
  /*  print out q-lms resolved PES spectrum */
  fp = fopen("peslm.dat","w") ;
  for ( ist = 0 ; ist <  nnm1fst ; ist++ )
    for ( iq = 0 ; iq < 3 ; iq++ ) 
      for ( iso = 0 ; iso < npesorb ; iso++ ) {
        sum = 0. ;
        for ( igstdeg = 0 ; igstdeg < gstdeg ; igstdeg++ ) {
          dum = pesmatele[iq][iso][igstdeg][ist] ;
          sum += dum*dum ;
	}
/*        if ( sum > EPSPES) */
	  fprintf(fp,"%11.6lf %16.10lf %2dq %2dn\n",
		  fstenergy[ist]-gstenergy, sum/gstdeg, iq-1, iso ) ;
      }
  fclose(fp) ;

/*   BUILD INTERMED STATE BASIS    {mstbas} */
  readconfs( nshells, lsh, &nconfs, &nelectrons, &occ, &nmstates ) ;
  makestates( nshells, lsh, sorb1sh, nconfs, nelectrons, occ, nmstates,
	      &mstbas) ;
  nhamele = calcham( &hamsparse, nmstates, mstbas, nconfs, nshells, lsh, 
		     sorb1sh, occ, ksish, e2p, e3d, cfdmat, hmag, thetamag ) ;
  for ( i = 0 ; i < nconfs ; i++ )
    free( occ[i] ) ;
  free( occ ) ;

/* diagonalize ham. and find intermed. eigen-states */
  ham = ( double * ) calloc( nmstates * nmstates, sizeof( double ) ) ;
  menergy = ( double * ) malloc( nmstates * sizeof( double ) ) ;
  for ( i = 0 ; i < nmstates * nmstates ; i++ ) ham[i] = 0. ;

  hamsparsep = hamsparse ; 
  for ( ist = 0 ; ist < nmstates ; ist++ ) {
    for ( k = 0 ; k < hamsparsep -> n ; k++ ) 
      ham[ ist + nmstates * hamsparsep -> j[k] ] = hamsparsep -> v[k] ;
    hamsparsep++ ;
  }
  spamadelete( hamsparse, nmstates ) ; 

  info =  solve_dsyev( nmstates, ham, menergy ) ; 
  if ( info ) { printf(" DSYEV : INFO = %d\n", info ) ; exit(1) ; } 

  nlevels = printspec( nmstates, menergy, &fst0deg, &fst0energy ) ; 

  printf("nlevels = %d\n", nlevels ) ;
  for ( i = 0 ; i < nmstates ; i++ ) { 
    printf(" %5.3lf", menergy[i]) ; 
    if ( i%20==19 ) printf("\n"); 
  }
  printf("\n"); 
/* define gamst */
  printf("Gamma grid: %8.3lf ", gamma0 ) ;
  for ( i = 0 ; i < ngam; i++ )  
    printf("%8.3lf %8.3lf ", egam[i], gam[i] ) ;
  printf("\n");
  gamst = (double *) malloc( nmstates*sizeof(double) );
  for ( jst = 0 ; jst < nmstates ; jst++ ) {
    gamma = gamma0 ; 
    if ( ngam > 0 ) {
      if ( menergy[jst] > egam[ngam-1] ) { 
        gamma = gam[ngam-1] ;
      }
      else {
        for ( i = 0 ; i < ngam-1 ; i++ ) 
	  if  ( egam[i] < menergy[jst] && menergy[jst] < egam[i+1] ) 
	    gamma = gam[i] ;
      }
    }
    gamst[jst] = gamma ;
  }
/*  for ( jst = 0 ; jst < nmstates ; jst++ ) 
    printf("%8.3lf %8.3lf\n", menergy[jst], gamst[jst] ) ;
*/
/* store intermed. eigen state vectors */
  mstvec = ( double ** ) malloc( nmstates * sizeof( double * ) ) ;
  for ( jst = 0 ; jst < nmstates ; jst++ )
    mstvec[jst] = ham + nmstates * jst ;    


/* DON'T print out intermediate eigensates 
  for ( jst = 0 ; jst < nmstates ; jst++ ) {
    printf("\n%8.5lf\n", menergy[jst]) ;  
    for ( i = 0 ; i < nmstates ; i++ ) { 
      if ( fabs( mstvec[jst][i] ) > 0.001 ) {
        printf("%8.5lf ", mstvec[jst][i]) ;
        fockdisplay( mstbas + i ) ;
        printf("\n") ;
      }
    }
  }
*/
/* construct 3d-t2g occupation operators for intermediate states */
/*
  po1plistitem = 0 ;
  nlistele = t2gocc2list( nshells, lsh, sorb1sh, 1, &po1plistitem ) ;
  printf("t2gocc Op. nlistele = %d.\n", nlistele ) ;
  t2gop = o1pmake( po1plistitem ) ;
  o1plistdelete( po1plistitem ) ;
  nhamele = o1ptospama( nmstates, mstbas, nmstates, mstbas,
			t2gop, &pt2gsmline0 ) ;
  po1plistitem = 0 ;
  nlistele = egocc2list( nshells, lsh, sorb1sh, 1, &po1plistitem ) ;
  printf("egocc Op. nlistele = %d.\n", nlistele ) ;
  egop = o1pmake( po1plistitem ) ;
  o1plistdelete( po1plistitem ) ;
  nhamele = o1ptospama( nmstates, mstbas, nmstates, mstbas,
			egop, &pegsmline0 ) ;
  po1plistitem = 0 ;
  nlistele = dxyocc2list( nshells, lsh, sorb1sh, 1, &po1plistitem ) ;
  printf("dxyocc Op. nlistele = %d.\n", nlistele ) ;
  dxyop = o1pmake( po1plistitem ) ;
  o1plistdelete( po1plistitem ) ;
  nhamele = o1ptospama( nmstates, mstbas, nmstates, mstbas,
			dxyop, &pdxysmline0 ) ;
*/
  /*
  for ( i = 0 ;  i < nmstates ; i++ ) {
    printf("\n<%3d|n( eg)|.> :", i ) ; spamalinewrite( pegsmline0[i] ) ;
    printf("       t2g     :" ) ;     spamalinewrite( pt2gsmline0[i] ) ;
  }
  */
  /*
  nt2g = ( double * ) malloc( nmstates * sizeof( double ) ) ;
  for ( jst = 0 ; jst < nmstates ; jst++ )
    nt2g[jst] = spamarealmatele( nmstates, mstvec[jst], nmstates, mstvec[jst], 
  			 pt2gsmline0 ) ;
  neg  = ( double * ) malloc( nmstates * sizeof( double ) ) ;
  for ( jst = 0 ; jst < nmstates ; jst++ )
    neg[jst] = spamarealmatele( nmstates, mstvec[jst], nmstates, mstvec[jst], 
  			 pegsmline0 ) ;
  ndxy = ( double * ) malloc( nmstates * sizeof( double ) ) ;
  for ( jst = 0 ; jst < nmstates ; jst++ )
    ndxy[jst] = spamarealmatele( nmstates, mstvec[jst], nmstates, mstvec[jst], 
  			 pdxysmline0 ) ;
  */

/* construct dipole operator  P^(1)_q = r C^(1)_q , q=-1,0,1 [Cowan (14.24)] */
  inishell = 0 ;
  for ( iq = 0 ; iq < 3 ; iq++ ) {
    q = iq - 1 ;
    po1plistitem = 0 ;
    nlistele = dipole2list( nshells, lsh, sorb1sh, inishell, q, radipmatele,
			       &po1plistitem) ;
    printf("Dipole Op. q = %d. nlistele = %d.\n", q, nlistele ) ;
    dipop[iq] = o1pmake( po1plistitem ) ;
    o1plistdelete( po1plistitem ) ;
    nhamele = o1ptospama( ngstbasis, gstbasis, nmstates, mstbas, 
			  dipop[iq], &pdipsmline0[iq] ) ;
    printf("numtotele = %d\n", nhamele ) ;
    /*    
      for ( i = 0 ;  i < nmstates ; i++ ) {
      printf("<%3d|r_{%2d}| . > . :", i, q ) ;   
      spamalinewrite( pdipsmline0[iq][i] ) ;
      }
    */
  }

  xasmatele = calloc3double( 3, gstdeg, nmstates ) ;
  for ( jst = 0 ; jst < nmstates ; jst++ )
    for ( iq = 0 ; iq < 3 ; iq++ )
      for ( igstdeg = 0 ; igstdeg < gstdeg ; igstdeg++ )
        xasmatele[iq][igstdeg][jst] =
        spamarealmatele( ngstbasis, gstvec[igstdeg], nmstates, mstvec[jst], 
  			 pdipsmline0[iq] ) ;

/*  print out q-resolved XAS spectrum */
  fp = fopen("xaq.dat","w") ;
  for ( i = 0 ; i <  nmstates ; i++ ) {
    sum = 0. ;
    for ( iq = 0 ; iq < 3 ; iq++ ) {
      dum = 0. ;
      for ( igstdeg = 0 ; igstdeg < gstdeg ; igstdeg++ )
        dum += xasmatele[iq][igstdeg][i] * xasmatele[iq][igstdeg][i] ;
      buf[iq] = dum / gstdeg ;
      sum += buf[iq] ;
    }
    if ( sum > EPSPES) {
      fprintf(fp,"%11.6lf %16.10lf", menergy[i]-gstenergy, sum ) ;
      for ( iq = 0 ; iq < 3 ; iq++ ) fprintf(fp," %16.10lf", buf[iq] ) ;
      fprintf(fp,"\n") ;
    }
  }   
  fclose(fp) ;

  /*  print out XAS for LCP RCP from positive x-axix (i.e. y(pi/2) rotation 
      LCP = [(-1)+sqrt2*(0)+(1)]/2 RCP = [(-1)-sqrt2*(0)+(1)]/2 */
  fp = fopen("xaqx.dat","w") ;
  for ( i = 0 ; i <  nmstates ; i++ ) {
    sum = 0. ; /* lcp */
    for ( igstdeg = 0 ; igstdeg < gstdeg ; igstdeg++ ) {
      dum = 0.5*( xasmatele[0][igstdeg][i] + sqrt(2)*xasmatele[1][igstdeg][i] + xasmatele[2][igstdeg][i] ) ;
      sum += dum*dum ;
    }
    buf[0] = sum/gstdeg ;
    sum = 0. ; /* lin.pol along x */
    for ( igstdeg = 0 ; igstdeg < gstdeg ; igstdeg++ ) {
      dum = sqrt(0.5)*( -xasmatele[0][igstdeg][i] + xasmatele[2][igstdeg][i] ) ;
      sum += dum*dum ;
    }
    buf[1] = sum/gstdeg ;
    sum = 0. ; /* rcp */
    for ( igstdeg = 0 ; igstdeg < gstdeg ; igstdeg++ ) {
      dum = 0.5*( xasmatele[0][igstdeg][i] - sqrt(2)*xasmatele[1][igstdeg][i] + xasmatele[2][igstdeg][i] ) ;
      sum += dum*dum ;
    }
    buf[2] = sum/gstdeg ;
    sum = buf[0]+buf[1]+buf[2] ;
    if ( sum > EPSPES) {
      fprintf(fp,"%11.6lf %16.10lf", menergy[i]-gstenergy, sum ) ;
      for ( iq = 0 ; iq < 3 ; iq++ ) fprintf(fp," %16.10lf", buf[iq] ) ;
      fprintf(fp,"\n") ;
    }
  }
  fclose(fp) ;

/*  print out t2g/eg/dxy-weighted XAS spectrum 
  fp = fopen("xte.dat","w") ;
  for ( i = 0 ; i <  nmstates ; i++ ) {
    sum = 0. ;
    for ( iq = 0 ; iq < 3 ; iq++ ) {
      dum = 0. ;
      for ( igstdeg = 0 ; igstdeg < gstdeg ; igstdeg++ )
        dum += xasmatele[iq][igstdeg][i] * xasmatele[iq][igstdeg][i] ;
      buf[iq] = dum / gstdeg ;
      sum += buf[iq] ;
    }
    if ( sum > EPSPES) {
      fprintf(fp,"%11.6lf %16.10lf %16.10lf %16.10lf %16.10lf", 
	menergy[i]-gstenergy, sum, sum*nt2g[i], sum*neg[i], sum*ndxy[i] ) ;
      for ( iq = 0 ; iq < 3 ; iq++ ) 
	fprintf(fp," %16.10lf %16.10lf", buf[iq]*nt2g[i], buf[iq]*neg[i] ) ;
     
      fprintf(fp,"\n") ;
    }
  }   
  fclose(fp) ;
*/

/* construct CVV AI operator  1=core, 2=valence, p=3rd-last shell
   VAI = Sum_{1,p,2!=2'}<1p|V|22'>c1+ cp+ c2 c2' */ 
  po2plistitem = 0 ;
  nlistele = vai2list( nshells, lsh, sorb1sh, &po2plistitem) ;
  printf("VAI. nlistele = %d.\n", nlistele ) ;
  vaiop = o2pmake( po2plistitem ) ;
  o2plistdelete( po2plistitem ) ;
  o2pprint( vaiop ) ;

/* calculate VAI matrix elements between basis states <Fbas|V|Ibas> */
  nhamele = o2ptospama( nmstates, mstbas, nfstates, fstate, 
	                vaiop, &pvaismline0 ) ;
  printf("numtotele = %d\n", nhamele ) ;
/*  for ( i = 0 ;  i < nfstates ; i++ ) {
    printf("<%3d|V| . > . :", i ) ; spamalinewrite( pvaismline0[i] ) ;  } 
*/

/* calculate VAI matrix elements <fist|V|M)=<i,iso|V|M) */
/* allocate fbasvmst */ /* here full matrix. could be improved using spama */
  fbasvmst = calloc2double( nfstates, nmstates ) ;
  for ( i = 0 ; i < nnm1fst ; i++ ) 
    for ( iso = 0 ; iso < npesorb ; iso++ ) {
      fist = iso*nnm1fst + i ;
      psml = pvaismline0 + fist ;
      for ( m = 0 ; m < nmstates ; m++ ) { 
        sum = 0. ;
        for ( k = 0 ; k < psml -> n ; k++ ) {
          if ( psml -> j[k] >= nmstates ) { printf("PANIC\n");exit(1);}
          sum  +=  psml -> v[k]  *  mstvec[m][ psml -> j[k] ] ;
        }
        fbasvmst[fist][m] = sum ;
      }
    }

/* calculate (F|V|M) as full matrix */
  k = 0 ;
  fstvmst = calloc2double( nfstates, nmstates ) ;
  for ( iso = 0 ; iso < npesorb ; iso++ ) 
    for ( j = 0 ; j < nnm1fst ; j++ ) {
      fjst = iso*nnm1fst + j ;
      for ( m = 0 ; m < nmstates ; m++ ) { 
        sum = 0. ;
        for ( i = 0 ; i < nnm1fst ; i++ ) 
	  sum += nm1fstvec[j][i] * fbasvmst[ iso*nnm1fst + i ][m] ; 
        fstvmst[fjst][m] = sum ;
        if ( fabs(sum) > EPSPES ) k++ ; 
      } 
    }  
  free2double( nfstates, fbasvmst ) ;       
  printf("Total # (F|V|M) = %d. #F=%d #M=%d NonZero/#F*#M=%lf\n", 
         k, nfstates, nmstates, k/(double)(nfstates*nmstates) ) ;

  rpesmatele = calloc4cplx( 3, npesorb, gstdeg, nnm1fst ) ;
  cmdum = ( double complex * ) malloc( nmstates * sizeof( double complex ) ) ;

  fp = fopen("rpes.dat","w") ;
  fpa = fopen("rpesalms.dat","w") ;
  fprintf( fpa,"%d %d %d %d\n", nomega, gstdeg, nnm1fst, npesorb ) ;
  fpy = fopen("rp.dat","w") ;
  fpc = fopen("rpc.dat","w") ;
  fpz = fopen("xmat.dat","w") ;
  fprintf( fpz,"%d %d x %d\n", nomega, npesorb/2, npesorb/2) ;
  
  fprintf(fpy,"%d\n", nnm1fst ) ;
  for ( j = 0 ; j < nnm1fst ; j++ ) 
    fprintf(fpy,"%11.6lf\n", gstenergy - nm1fenergy[j] ) ;
  fprintf(fpy,"%d\n", nomega ) ;

  fprintf(fpc,"%d\n", nnm1fst ) ;
  for ( j = 0 ; j < nnm1fst ; j++ ) 
    fprintf(fpc,"%11.6lf\n", gstenergy - nm1fenergy[j] ) ;
  fprintf(fpc,"%d\n", nomega ) ;

  fpx = fopen("xaqc.dat","w") ;
/* main calculation: loop over omega */
  for ( omega = ommin ; omega <= ommax+EPSILON ; omega += deltaom ) {
/*    fprintf(fp,"omega=%lf\n",omega); */

/* calculate resonant matrix element (F|V|M)(M|D|G)/[om+EG-EM+iGa] */
  for ( iq = 0 ; iq < 3 ; iq++ ) {
    csumxaq[iq] = 0. ; 
    for ( igstdeg = 0 ; igstdeg < gstdeg ; igstdeg++ ) {
      for ( m = 0 ; m < nmstates ; m++ )  {
	/* linear increase: 
	gamst[m] = 0.2 + ((menergy[m]>-5.)?(menergy[m]+5.):0.)*0.04 ; 
        */
        cmdum[m] = xasmatele[iq][igstdeg][m]
                 / ( omega + gstenergy - menergy[m] + I*gamst[m] ) ;
        csumxaq[iq] += xasmatele[iq][igstdeg][m]*cmdum[m] ;
      } 
      for ( iso = 0 ; iso < npesorb ; iso++ ) {
        for ( j = 0 ; j < nnm1fst ; j++ ) {
          fjst = iso*nnm1fst + j ;
          csum = 0. ;
          for ( m = 0 ; m < nmstates ; m++ ) 
            csum += fstvmst[fjst][m]*cmdum[m] ;
          rpesmatele[iq][iso][igstdeg][j] = csum
                 + pesmatele[iq][iso][igstdeg][j] ;
        }
      }
    }
  }
  printrpesmatele( fpa, omega, npesorb, nnm1fst, gstdeg, gstweight,
		   nm1fenergy, gstenergy, rpesmatele ) ;
  printrpesxmat( fpz, omega, npesorb, nnm1fst, gstdeg, gstweight, rpesmatele );

  
/*  print out total RPES spectrum Sum_lms LCP and RCP */
  fprintf(fpc,"\n%15.8e\n\n", omega) ;
  for ( j = 0 ; j < nnm1fst ; j++ ) {
    suml = sum = sumr = 0. ;
    for ( iso = 0 ; iso < npesorb ; iso++ ) {
      for ( igstdeg = 0 ; igstdeg < gstdeg ; igstdeg++ ) {
	for ( zdum = 0., iq = 0 ; iq < 3 ; iq++ )
	  zdum += eps[0][iq] * rpesmatele[iq][iso][igstdeg][j] ;
	suml += conj(zdum)*zdum ;
	for ( zdum = 0., iq = 0 ; iq < 3 ; iq++ )
          zdum += eps[1][iq] * rpesmatele[iq][iso][igstdeg][j] ;
        sum  += conj(zdum)*zdum ;
	for ( zdum = 0., iq = 0 ; iq < 3 ; iq++ )
	  zdum += eps[2][iq] * rpesmatele[iq][iso][igstdeg][j] ;
	sumr += conj(zdum)*zdum ;
      }
    }
    fprintf(fpc,"%15.8e %15.8e %15.8e\n",suml/gstdeg,sum/gstdeg,sumr/gstdeg);
  }
/*  print out total RPES spectrum Sum_q Sum_lms */
  fprintf(fpy,"\n%15.8e\n\n", omega) ;
  for ( j = 0 ; j < nnm1fst ; j++ ) {
    sum = 0. ;
    for ( iq = 0 ; iq < 3 ; iq++ ) {      
      for ( iso = 0 ; iso < npesorb ; iso++ ) {
        for ( igstdeg = 0 ; igstdeg < gstdeg ; igstdeg++ ) {
          zdum = rpesmatele[iq][iso][igstdeg][j] ;
          sum += conj(zdum)*zdum ; 
        }
      }
    }
    if ( sum > EPSPES) 
      fprintf(fp,"%10.6lf %11.6lf %15.8e\n",
              omega, gstenergy - nm1fenergy[j], sum/gstdeg ) ;
    fprintf(fpy,"%15.8e\n", sum/gstdeg ) ;    
  }
  /* print broadened XAS spectrum */
  fprintf(fpx,"%10.6lf ", omega ) ;
  sum = 0. ;
  for ( iq = 0 ; iq < 3 ; iq++ ) {
    dum = -cimag( csumxaq[iq] ) / gstdeg / PI ;
    fprintf(fpx,"%15.8e ", dum ) ; 
    sum += dum ;
  }
  fprintf(fpx,"%15.8e\n", sum) ;

  } /* end omega loop */
  fclose(fpz) ;
  fclose(fpy) ;
  fclose(fpa) ;
  fclose(fpx) ;
  fclose(fp) ;
  free( cmdum ) ;  
  free4cplx( 3, npesorb, gstdeg, rpesmatele ) ;



  for ( iq = 0 ; iq < 3 ; iq++ ) 
    spamadelete( pdipsmline0[iq], nmstates ) ;
  for ( igstdeg = 0 ; igstdeg < gstdeg ; igstdeg++ )
    free( gstvec[igstdeg] ) ;
  free( gstvec ) ;
  free( gstbasis ) ;
}


/* definition of functions at top */

void printstatevector( int nstates, double *vec, struct Fock *state ) 
{
  int i ;
  for ( i = 0 ; i < nstates ; i++ ) {
    /* long version */
    if ( fabs(vec[i]) > EPSILON ) {
      printf("%15.10lf ", vec[i] ) ;
      fockdisplay( &state[i] ) ;
      printf("\n") ;
    }
  }
  printf("\n") ;
}

int printspec( int dim, double *lambda, int *pgstdeg, double *pgstenergy ) 
{
    int nlevels, degeneracy, i, gstdeg ;
    double gstenergy, newlambda, oldlambda ;
    nlevels = 0 ;
    for ( i = 0 ; i < dim ; i++ ) {
      newlambda = lambda[i] ;
    if ( nlevels == 0 ) {
      oldlambda = newlambda ;
      degeneracy = 1 ;
      nlevels = 1 ;
    } 
    else if ( fabs( newlambda - oldlambda ) < EPSILON )
      degeneracy++ ;
    else {
      printf(" %12.6lf  -  degeneracy = %d\n", oldlambda, degeneracy ) ;
      if ( nlevels == 1 ) {
	gstenergy = oldlambda ;
	gstdeg = degeneracy ;
      }
      oldlambda = newlambda ;
      degeneracy = 1 ;
      nlevels++ ;
    }
  }      
  printf(" %12.6lf  -  degeneracy = %d\n", oldlambda, degeneracy ) ;
  if ( nlevels == 1 ) {
    gstenergy = oldlambda ;
    gstdeg = degeneracy ;
  }

  *pgstdeg = gstdeg ;
  *pgstenergy = gstenergy ;
  return nlevels ;
}
void printrpesmatele( FILE *fp, double omega, int npesorb, int nnm1fst,
		      int gstdeg, double *gstweight, double *nm1fenergy,
		      double gstenergy, complex ****rpesmatele )
{
  int j, iso, iq, igstdeg ;
  complex dum ;
  fprintf(fp,"%17.10e\n", omega ) ;
  for ( igstdeg = 0 ; igstdeg < gstdeg ; igstdeg++ ) {
    fprintf(fp,"%lf\n", gstweight[igstdeg] ) ;	    
    /* fprintf(fp,"# ground state No. %d / %d\n", igstdeg+1, gstdeg ) ;	*/    
    for ( j = 0 ; j < nnm1fst ; j++ ) {
      fprintf(fp,"%17.10e\n", gstenergy - nm1fenergy[j] ) ;
      for ( iso = 0 ; iso < npesorb ; iso++ ) {
	for ( iq = 0 ; iq < 3 ; iq++ ) {
	  dum = rpesmatele[iq][iso][igstdeg][j] ;
	  fprintf(fp,"   %13.6e %13.6e", creal(dum), cimag(dum) ) ;
	}
	fprintf(fp,"\n") ;	
      }
    }
  }
}  
void printrpesxmat( FILE *fp, double omega, int npesorb, int nnm1fst,
		    int gstdeg, double *gstweight, complex ****rpesmatele )
{
  int io, jo, iq, igstdeg, j, isp, npeorb ;
  complex sum, sumg ;
  npeorb = npesorb/2 ;
  fprintf(fp,"%17.10e\n", omega ) ;
  for ( io = 0 ; io < npeorb ; io++ ) {
    for ( jo = 0 ; jo < npeorb ; jo++ ) {
      for ( iq = 0 ; iq < 3 ; iq++ ) {
	sum = 0. ;
	for ( igstdeg = 0 ; igstdeg < gstdeg ; igstdeg++ ) {
	  sumg = 0.;
	  for ( j = 0 ; j < nnm1fst ; j++ ) 
	    for ( isp = 0 ; isp < 2 ; isp++ ) 
	      sumg  +=  rpesmatele[iq][2*io+isp][igstdeg][j]
		* conj( rpesmatele[iq][2*jo+isp][igstdeg][j] ) ;
	  sum += gstweight[igstdeg] * sumg ;
	}
	fprintf(fp," %13.6e %13.6e", creal(sum), cimag(sum) ) ;
      }
      fprintf(fp,"\n") ;	
    }
  }
  fprintf(fp,"\n") ;	
}

/*  print out angle-indep. RPES spectrum
  for ( iq = 0 ; iq < 3 ; iq++ ) {  
   for ( iqp = iq ; iqp < 3 ; iqp++ ) {
 
    for ( j = 0 ; j < nnm1fst ; j++ ) {
     sum = 0. ;
      for ( iso = 0 ; iso < npesorb ; iso++ ) {
        for ( igstdeg = 0 ; igstdeg < gstdeg ; igstdeg++ ) {
          zdum = rpesmatele[iq][iso][igstdeg][j] ;
          sum += conj(zdum)*zdum ; 
        }
      }
    }
    if ( sum > EPSPES) 
      fprintf(fp,"%10.6lf %11.6lf %15.8lf\n",
              omega, nm1fenergy[j]-gstenergy, sum/gstdeg ) ;
  }
*/
