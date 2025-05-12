#include <stdio.h>
#include <stdlib.h>
#include "globals.h"

void readshells( int *pnshells, int **plsh, int **psorb1sh, double **pksish) 
{
  int nshells, *lsh, *sorb1sh, i ;
  double *ksish ;

  printf("nshells = ") ; scanf("%d", &nshells ) ;
  lsh     = ( int * ) malloc( nshells * sizeof( int ) ) ;
  sorb1sh = ( int * ) malloc( ( nshells + 1 ) * sizeof( int ) ) ;
  ksish = ( double * ) malloc( nshells * sizeof( double ) ) ;

  sorb1sh[0] = 0 ;
  printf("Enter l-values for all %d shells: ", nshells ) ;
  for ( i = 0 ; i < nshells ; i++ ) {
    scanf("%d", &lsh[i] ) ;
    sorb1sh[i+1] = sorb1sh[i] + 4 * lsh[i] + 2 ;
  }
  printf("Enter ksi-values for all %d shells: ", nshells ) ;
  for ( i = 0 ; i < nshells ; i++ ) 
    scanf("%lf", &ksish[i] ) ;
  
  printf("nshells = %d, nsorbs = %d\n",	nshells, sorb1sh[nshells] ) ;
  for ( i = 0 ; i < nshells ; i++ ) 
    printf("%3d %3d\n", lsh[i], sorb1sh[i+1] ) ; 

  *pnshells = nshells ;
  *plsh = lsh ;
  *psorb1sh = sorb1sh ;
  *pksish = ksish ;
}

void readconfs( int nshells, int *lsh, int *pnconfs, int *pnelectrons,
		int ***pocc, int *pnstates ) 
{
  int nconfs, nelectrons, nelectronsbefore, nstatesinconf, i, j;
  int **occ, nstates ;

  printf("nconfs = ") ; scanf("%d", &nconfs ) ;
  occ = ( int ** ) malloc( nconfs * sizeof( int * ) ) ;

  nstates = 0 ;
  for ( i = 0 ; i < nconfs ; i++ ) {
    occ[i] = ( int * ) malloc( nshells * sizeof( int ) ) ;

    printf("Conf No %d. Enter occupancy for all %d shells: ", i, nshells ) ;

    nelectrons = 0 ;
    nstatesinconf = 1 ;
    for ( j = 0 ; j < nshells ; j++ ) {
      scanf("%d", &occ[i][j] ) ;
      if ( occ[i][j] < 0 || occ[i][j] > 4 * lsh[j] + 2 ) {
	printf("occ[i][j] < 0 || occ[i][j] > 4 * lsh[j] + 2\n" ) ; exit(1) ; 
      }
      nelectrons += occ[i][j] ;
      nstatesinconf *= noverk( 4 * lsh[j] + 2, occ[i][j] ) ;
    }
    if ( i > 0 && nelectrons != nelectronsbefore ) {
      printf("\n\n nelectrons != nelectronsbefore\n\n") ; exit(1) ; 
    }
    nelectronsbefore = nelectrons ;
    nstates +=  nstatesinconf ;
  }
  *pnconfs = nconfs ;
  *pnelectrons = nelectrons ;
  *pocc = occ ;
  *pnstates = nstates ;
}
