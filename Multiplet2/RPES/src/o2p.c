#include <stdio.h>
#include <stdlib.h>
#include "globals.h"

void o2plistitemadd( struct O2plistitem **ppitem0, 
		      int i1, int i2, int i3, int i4, double v )
{
  struct O2plistitem *pitem ;
  pitem = ( struct O2plistitem * ) malloc( sizeof( struct O2plistitem ) ) ;
  pitem -> i1 = i1 ;
  pitem -> i2 = i2 ;
  pitem -> i3 = i3 ;
  pitem -> i4 = i4 ;
  pitem -> v  = v ;
  pitem -> next = *ppitem0 ;
  *ppitem0 = pitem ;
}


void o2plistdelete( struct O2plistitem *pitem )
{
  int counter = 0 ;
  struct O2plistitem *pnextitem ;
  /*  printf("Deleting o2plist. Deleted items :\n") ; */
  while ( pitem != 0 ) {
    /*
    printf("(%d %d %d %d %lf) ", pitem -> i1, pitem -> i2, 
	   pitem -> i3, pitem -> i4, pitem -> v ) ;
    */
    pnextitem = pitem -> next ;
    free( pitem ) ;
    pitem = pnextitem ;
    counter++ ;
  }
  /*  printf("\n In all %d items have been deleted.\n", counter ) ; */
}


int o2pnfindi1( struct O2plistitem *listp, int itab[], int ntab )
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


int o2pnfindi2( struct O2plistitem *listp, int itab[], int ntab, int i1)
{
  int i, n ;

  for ( i = 0 ; i < ntab ; i++ ) 
    itab[i] = 0 ;

  while ( listp != 0 ) {
    if ( listp->i1 == i1 ) { 
      itab[ listp->i2 ] = 1 ;
    }
    listp = listp -> next ;
  } 		

  n = 0 ;
  for ( i = 0 ; i < ntab ; i++ ) 
    if ( itab[i] == 1 ) 
      n++ ;
 
  return n ;
}  


int o2pnfindi3( struct O2plistitem *listp, int itab[], int ntab,
		 int i1, int i2 )
{
  int i, n ;

  for ( i = 0 ; i < ntab ; i++ ) 
    itab[i] = 0 ;

  while ( listp != 0 ) {
    if ( listp->i1 == i1 && listp->i2 == i2 ) { 
      itab[ listp->i3 ] = 1 ;
    }
    listp = listp -> next ;
  } 		

  n = 0 ;
  for ( i = 0 ; i < ntab ; i++ ) 
    if ( itab[i] == 1 ) 
      n++ ;
 
  return n ;
}  


int o2pnfindi4v( struct O2plistitem *listp, int itab[], double vtab[],
	          int ntab, int i1, int i2, int i3 )
{
  int i, n ;

  for ( i = 0 ; i < ntab ; i++ ) {
    itab[i] = 0  ;
    vtab[i] = 0. ;
  }

  while ( listp != 0 ) {
    if ( listp->i1 == i1 && listp->i2 == i2 && listp->i3 == i3 ) { 
      itab[ listp->i4 ] = 1 ;
      vtab[ listp->i4 ] += listp->v ;
    }
    listp = listp -> next ;
  } 		

  n = 0 ;
  for ( i = 0 ; i < ntab ; i++ ) 
    if ( itab[i] == 1 ) 
      n++ ;
 
  return n ;
}  

struct O2p o2pmake( struct O2plistitem *listp0 )
{
  int itab[NBITS], i, k, k1, k2, k3,  n1, n2, n3, n4 ;
  double vin, vtab[NBITS] ;
  struct O2p op2p ;

  n1 = o2pnfindi1( listp0, itab, NBITS) ;
  op2p.n = n1 ;
  op2p.i = ( int * ) malloc( n1 * sizeof( int ) ) ; 
  op2p.p = ( struct O2p2 * ) malloc( n1 * sizeof( struct O2p2 ) ) ;
  for ( k = 0, i = 0 ; i < NBITS ; i++ )
    if ( itab[i] == 1 ) 
      op2p.i[k++] = i ;
  for ( k1 = 0 ; k1 < n1 ; k1++ ) {
    n2 = o2pnfindi2( listp0, itab, NBITS, op2p.i[k1] ) ;
    op2p.p[k1].n = n2 ;
    op2p.p[k1].i = ( int * ) malloc( n2 * sizeof( int ) ) ;
    op2p.p[k1].p = ( struct O2p3 * ) malloc( n2 * sizeof( struct O2p3 ) ) ;
    for ( k = 0, i = 0 ; i < NBITS ; i++ )
      if ( itab[i] == 1 ) {
	op2p.p[k1].i[k] = i ;
	k++ ;
      }
    for ( k2 = 0 ; k2 < n2 ; k2++ ) {
      n3 = o2pnfindi3( listp0, itab, NBITS, op2p.i[k1], op2p.p[k1].i[k2] ) ;
      op2p.p[k1].p[k2].n = n3 ;
      op2p.p[k1].p[k2].i = ( int * ) malloc( n3 * sizeof( int ) ) ;
      op2p.p[k1].p[k2].p =
	( struct O2p4 * ) malloc( n3 * sizeof( struct O2p4 ) ) ;
      for ( k = 0, i = 0 ; i < NBITS ; i++ )
	if ( itab[i] == 1 ) {
	  op2p.p[k1].p[k2].i[k] = i ;
	  k++ ;
	}
      for ( k3 = 0 ; k3 < n3 ; k3++ ) {
	n4 = o2pnfindi4v( listp0, itab, vtab, NBITS, op2p.i[k1], 
			   op2p.p[k1].i[k2], op2p.p[k1].p[k2].i[k3] ) ;
	op2p.p[k1].p[k2].p[k3].n = n4 ;
	op2p.p[k1].p[k2].p[k3].i = ( int * ) malloc( n4 * sizeof( int ) ) ;
	op2p.p[k1].p[k2].p[k3].v = ( double *) malloc( n4 * sizeof( double ) );
	for ( k = 0, i = 0 ; i < NBITS ; i++ )
	  if ( itab[i] == 1 ) {
	    op2p.p[k1].p[k2].p[k3].i[k] = i ;
	    op2p.p[k1].p[k2].p[k3].v[k] = vtab[i] ;
	    k++ ;
	  }
      }
    }
  }
  return op2p ;
}


void o2pprint( struct O2p vcb ) {
  int k1, k2, k3, k4 ;
/* print out the dimensions n1, n2(i1), n3(i1,i2), n4(i1,i2,i3,i4)  */
/*
  printf("Dimensions. n1 = %d\n",  vcb.n ) ;
  for ( k1 = 0 ; k1 < vcb.n ; k1++ ) {
    printf("i1 = %d, n2 = %d \n",  vcb.i[k1], vcb.p[k1].n ) ;
    for ( k2 = 0 ; k2 < vcb.p[k1].n ; k2++ ) {
      printf("i2 = %d, n3 = %d \n",  vcb.p[k1].i[k2], vcb.p[k1].p[k2].n ) ;
      for ( k3 = 0 ; k3 < vcb.p[k1].p[k2].n ; k3++ ) 
	printf("i3 = %d, n4 = %d ; ",  
	       vcb.p[k1].p[k2].i[k3],  vcb.p[k1].p[k2].p[k3].n ) ;
      printf("\n") ;
    }
  }
*/
/* print out the dynamical 4-dim array v2p */
  for ( k1 = 0 ; k1 < vcb.n ; k1++ ) 
    for ( k2 = 0 ; k2 < vcb.p[k1].n ; k2++ )
      for ( k3 = 0 ; k3 < vcb.p[k1].p[k2].n ; k3++ )
	for ( k4 = 0 ; k4 < vcb.p[k1].p[k2].p[k3].n ; k4++ )
	  printf("%2d %2d %2d %2d %11.6lf\n",  vcb.i[k1], vcb.p[k1].i[k2],
		 vcb.p[k1].p[k2].i[k3], vcb.p[k1].p[k2].p[k3].i[k4],
		 vcb.p[k1].p[k2].p[k3].v[k4] ) ;

}

void o2pdelete( struct O2p vcb ) {
  int k1, k2, k3 ;
  for ( k1 = 0 ; k1 < vcb.n ; k1++ ) {
    for ( k2 = 0 ; k2 < vcb.p[k1].n ; k2++ ) {
      for ( k3 = 0 ; k3 < vcb.p[k1].p[k2].n ; k3++ ) {
	free( vcb.p[k1].p[k2].p[k3].i ) ;
	free( vcb.p[k1].p[k2].p[k3].v ) ;
      }
      free( vcb.p[k1].p[k2].i ) ;
      free( vcb.p[k1].p[k2].p ) ;
    }
    free( vcb.p[k1].i ) ;
    free( vcb.p[k1].p ) ;
  }
  free( vcb.i ) ;
  free( vcb.p ) ;
}


int o2ptospama( int nketbasis, struct Fock *pketbasis, int nbrabasis, 
		struct Fock *pbrabasis, struct O2p op2p, 
		struct Spamaline **ppsml0 ) 
{
  int ketj, brai, i1, i2, i3, i4, k1, k2, k3, k4, sign1, sign2, sign3, sign4 ;
  int totnumele ;
  double *columnj, val ;
  struct Fock state0, state1, state2, state3, state4 ;
  struct O2p2 op2p2 ;
  struct O2p3 op2p3 ;
  struct O2p4 op2p4 ;
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

/* apply op2p to state */
    for ( k1 = 0 ; k1 < op2p.n ; k1++ ) {
      state1 = state0 ;
      i1 = op2p.i[k1] ;
      sign1 = fockdes( &state1, i1 ) ;
      if ( sign1 != 0 ) {
	op2p2 = op2p.p[k1] ;
	for ( k2 = 0 ; k2 < op2p2.n ; k2++ ) {
	  state2 = state1 ;
	  i2 = op2p2.i[k2] ;
	  sign2 = fockdes( &state2, i2 ) ;
	  if ( sign2 != 0 ) {
	    op2p3 = op2p.p[k1].p[k2] ;
	    for ( k3 = 0 ; k3 < op2p3.n ; k3++ ) {
	      state3 = state2 ;
	      i3 = op2p3.i[k3] ;
	      sign3 = fockcre( &state3, i3 ) ;
	      if ( sign3 != 0 ) {
		op2p4 = op2p.p[k1].p[k2].p[k3] ;
		for ( k4 = 0 ; k4 < op2p4.n ; k4++ ) {
		  state4 = state3 ;
		  i4 = op2p4.i[k4] ;
		  sign4 = fockcre( &state4, i4 ) ;
		  if ( sign4 != 0 ) {
		    val = sign4 * sign3 * sign2 * sign1 *  op2p4.v[k4] ;
		    if( findstate( nbrabasis, pbrabasis, &state4, &brai ) ) {
		      columnj[brai] += val ;
		    }
		    else { 
		      printf("c+%d c+%d c%d c%d  -> %lf ", i4,i3,i2,i1, val );
		      fockdisplay( &state4 ) ; 
		      printf(" = new state under No %d\n", brai ) ;
		    }
		  }
		}
	      }
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
