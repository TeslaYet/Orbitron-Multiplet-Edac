#include <stdio.h>
#include <stdlib.h>
#include "globals.h"


void o1plistitemadd( struct O1plistitem **ppitem0, int i1, int i2, double v )
{
  struct O1plistitem *pitem ;
  pitem = ( struct O1plistitem * ) malloc( sizeof( struct O1plistitem ) ) ;
  pitem -> i1 = i1 ;
  pitem -> i2 = i2 ;
  pitem -> v  = v ;
  pitem -> next = *ppitem0 ;
  *ppitem0 = pitem ;
}


void o1plistdelete( struct O1plistitem *pitem )
{
  int counter = 0 ;
  struct O1plistitem *pnextitem ;
  printf("Deleting op1plist. Deleted items :\n") ;
  while ( pitem != 0 ) {
    printf("(%d,%d,%7lf) ", pitem -> i1, pitem -> i2, pitem -> v ) ;
    pnextitem = pitem -> next ;
    free( pitem ) ;
    pitem = pnextitem ;
    counter++ ;
  }
  printf("\n In all %d items have been deleted.\n", counter ) ;
}


int o1pnfindi1( struct O1plistitem *listp, int itab[], int ntab )
{
  int i, n ;

  /* clear itab */
  for ( i = 0 ; i < ntab ; i++ ) 
    itab[i] = 0 ;

  /* find all different i's */
  while ( listp != 0 ) {
    itab[ listp->i1 ] = 1 ;
    listp = listp -> next ;
  } 

  /* count i's */
  n = 0 ;
  for ( i = 0 ; i < ntab ; i++ ) 
    if ( itab[i] == 1 ) 
      n++ ;
 
  /* printout */  
/*
  printf("n = %d. i1 = ", n ) ;
  for ( i = 0 ; i < ntab ; i++ )
    if ( itab[i] == 1 ) 
      printf("%d ", i ) ;
  printf("\n") ;
*/

  return n ;
}  

int o1pnfindi2v( struct O1plistitem *listp, int itab[], double vtab[],
	          int ntab, int i1 )
{
  int i, n ;

  for ( i = 0 ; i < ntab ; i++ ) {
    itab[i] = 0  ;
    vtab[i] = 0. ;
  }

  while ( listp != 0 ) {
    if ( listp->i1 == i1 ) { 
      itab[ listp->i2 ] = 1 ;
      vtab[ listp->i2 ] += listp->v ;
    }
    listp = listp -> next ;
  } 		

  n = 0 ;
  for ( i = 0 ; i < ntab ; i++ ) 
    if ( itab[i] == 1 ) 
      n++ ;
 
  return n ;
}  

struct O1p o1pmake( struct O1plistitem *listp0 )
{
  int itab[NBITS], i, k, k1, n1, n2 ;
  double vin, vtab[NBITS] ;
  struct O1p op1p ;

  n1 = o1pnfindi1( listp0, itab, NBITS) ;
  op1p.n = n1 ;
  op1p.i = ( int * ) malloc( n1 * sizeof( int ) ) ; 
  op1p.p = ( struct O1p2 * ) malloc( n1 * sizeof( struct O1p2 ) ) ;
  for ( k = 0, i = 0 ; i < NBITS ; i++ )
    if ( itab[i] == 1 ) 
      op1p.i[k++] = i ;
  for ( k1 = 0 ; k1 < n1 ; k1++ ) {
    n2 = o1pnfindi2v( listp0, itab, vtab, NBITS, op1p.i[k1] ) ;
    op1p.p[k1].n = n2 ;
    op1p.p[k1].i = ( int * ) malloc( n2 * sizeof( int ) ) ;
    op1p.p[k1].v = ( double *) malloc( n2 * sizeof( double ) ) ;
    for ( k = 0, i = 0 ; i < NBITS ; i++ )
      if ( itab[i] == 1 ) {
	op1p.p[k1].i[k] = i ;
	op1p.p[k1].v[k] = vtab[i] ;
	k++ ;
      }
  }
  return op1p ;
}


void o1pprint( struct O1p h1p ) {
  int k1, k2 ;
/* print out the dynamical 2-dim array v1p */
  for ( k1 = 0 ; k1 < h1p.n ; k1++ ) 
    for ( k2 = 0 ; k2 < h1p.p[k1].n ; k2++ )
      printf("%d %d %lf\n",  h1p.i[k1], h1p.p[k1].i[k2],
	     h1p.p[k1].v[k2] ) ;
}


void o1pdelete( struct O1p h1p ) {
  int k1 ;
  for ( k1 = 0 ; k1 < h1p.n ; k1++ ) {
    free( h1p.p[k1].i ) ;
    free( h1p.p[k1].v ) ;
  }
  free( h1p.i ) ;
  free( h1p.p ) ;
}


int o1ptospama( int nketbasis, struct Fock *pketbasis, int nbrabasis, 
		struct Fock *pbrabasis, struct O1p op1p, 
		struct Spamaline **ppsml0 ) 
{
  int ketj, brai, i1, i2, k1, k2, sign1, sign2, totnumele ;
  double *columnj, val ;
  struct Fock state0, state1, state2 ;
  struct O1p2 op1p2 ;
  struct Spamaline *psmcol0 ;
	

  psmcol0 = ( struct Spamaline * )
    malloc( nketbasis * sizeof( struct Spamaline ) ) ;

  columnj = ( double * ) malloc( nbrabasis * sizeof( double ) ) ;
  totnumele = 0 ;

  for ( ketj = 0 ; ketj < nketbasis ; ketj++ ) {
    state0 = pketbasis[ketj] ;
    /* fockdisplay( &state0 ) ; printf("\n") ; */

    /* clear columnj */
    for ( brai = 0 ; brai < nbrabasis ; brai++ ) 
      columnj[brai] = 0. ;

/* apply op1p to state */
    for ( k1 = 0 ; k1 < op1p.n ; k1++ ) {
      state1 = state0 ;
      i1 = op1p.i[k1] ;
      sign1 = fockdes( &state1, i1 ) ;
      if ( sign1 != 0 ) {
	op1p2 = op1p.p[k1] ;
	for ( k2 = 0 ; k2 < op1p2.n ; k2++ ) {
	  state2 = state1 ;
	  i2 = op1p2.i[k2] ;
	  sign2 = fockcre( &state2, i2 ) ;
	  if ( sign2 != 0 ) {
	    val = sign2 * sign1 *  op1p2.v[k2] ;
	    /*
	       printf("c+%d c%d  -> %lf ", i2, i1, val ); 
	       fockdisplay( &state2 ) ;
	    */ 
	    if( findstate(nbrabasis, pbrabasis, &state2, &brai ) ) {
	      /* printf(" = state No %d\n", brai ) ;  */
	       columnj[brai] += val ;
	    }
	    else { 
	      printf("c+%d c%d  -> %lf ", i2, i1, val ); 
	      fockdisplay( &state2 ) ;
	      printf(" = new state under No %d\n", brai ) ;
	    } 
	  }
	}
      }
    }
    totnumele += darray2spamaline( psmcol0 + ketj, nbrabasis, columnj ) ;
  }


  free( columnj ) ; 
/*
  printf("\n psmcol0[] totnumele = %d.\n\n", totnumele ) ;
  
  for ( ketj = 0 ; ketj < nketbasis ; ketj++ ) {
    printf("ketj =%3d --", ketj ) ;   
    spamalinewrite( psmcol0[ketj] ) ;
  }
*/

  totnumele = spamatranspose( &psmcol0, nketbasis, nbrabasis ) ;

/*
  for ( brai = 0 ; brai < nbrabasis ; brai++ ) {
    printf("brai =%3d --", brai) ;   
    spamalinewrite( psmcol0[brai] ) ;
  }
*/

  *ppsml0 = psmcol0 ;
  return totnumele ;
}
