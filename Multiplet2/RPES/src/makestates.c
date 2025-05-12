#include <stdio.h>
#include <stdlib.h>
#include "globals.h"

void makestates( int nshells, int *lsh, int *sorb1sh, int nconfs, 
		int nelectrons, int **occ, int nstates0, 
		struct Fock **pstate) 
{
  int nstates, i, k, ish, iconf, ishrun, iel, jel, ist ;
  int *nstsh, ***combsh, *kst, *occorbs ;
  struct Fock *state ;

  printf("nconfs = %d, nstates = %d\n",	nconfs, nstates0 ) ;
  state = ( struct Fock * ) calloc( nstates0, sizeof( struct Fock ) ) ;

  combsh  = ( int *** ) malloc( nshells * sizeof( int ** ) ) ;
  nstsh = ( int * ) malloc( nshells * sizeof( int ) ) ;
  
  occorbs = ( int * ) malloc( nelectrons * sizeof( int ) ) ;
  kst = ( int * ) malloc( nshells * sizeof( int ) ) ;

  nstates = 0 ;
  for ( iconf = 0 ; iconf < nconfs ; iconf++ ) {
    printf("Conf No %2d = ( ", iconf ) ;
    for ( ish = 0 ; ish < nshells ; ish++ )
      printf("%d ", occ[iconf][ish] ) ;
    printf(").\n") ;
    for ( ish = 0 ; ish < nshells ; ish++ ) {
      k = occ[iconf][ish] ;
      nstsh[ish] = 
	combinations( 4*lsh[ish]+2, k , sorb1sh[ish], &combsh[ish] ) ;
    }
    for ( i = 0 ; i < nshells ; i++ )
      kst[i] = 0 ;   /* clear kst[] */
    do { 
      iel = 0 ;
      for ( ish = 0 ; ish < nshells ; ish++ )
	for ( jel = 0 ; jel < occ[iconf][ish] ; jel++ ) 
	  occorbs[iel++] = combsh[ish][kst[ish]][jel] ;
      ishrun = 0 ;
      while ( ishrun < nshells && kst[ishrun] == nstsh[ishrun] - 1 ) 
	ishrun++ ;
      if ( ishrun < nshells ) kst[ishrun]++ ;
      for ( ish = 0 ;  ish < ishrun ; ish++ ) 
	kst[ish] = 0 ;
      printf("%4d : ", nstates ) ;
      for ( jel = 0 ; jel < nelectrons ;  jel++ ) 
	printf("%2d ", occorbs[jel] ) ;
      printf("\n") ;
      fockinit( state + nstates, iel, occorbs ) ; 
      nstates++ ;
    } while ( ishrun < nshells ) ; 
    /* instead of the " do { ... } while " above, one may as well use: */
    /*    nstatesinconf = state1conf[iconf+1] - state1conf[iconf] ;    */
    /*    for ( ist = 0 ; ist < nstatesinconf ; ist++ ) {  ...      }  */

    /* deallocate (**) memory of this awful combsh[][*][*] */
    for ( ish = 0 ; ish < nshells ; ish++ ) {
      for ( i = 0 ; i < nstsh[ish] ; i++ )
	free( combsh[ish][i] ) ;
      free( combsh[ish] ) ;
    }
  }
  free( kst ) ;
  free( occorbs ) ;
  free( nstsh ) ;
  free( combsh ) ;


  if ( nstates != nstates0 ) { printf("nstates != nstates0\n") ; exit(1) ; }

  /* sorting the states */
  focklinearsort( nstates, state ) ;
  for ( ist = 0 ; ist < nstates ; ist++ ) {
    printf("%3d:%18llu)=", ist, state[ist].n) ; 
    fockdisplay( state + ist ) ;
    printf("\n") ;
  }
  *pstate = state ;
  /*  return nstates ;
   */
}
