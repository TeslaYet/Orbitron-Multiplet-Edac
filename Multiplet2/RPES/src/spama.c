#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "globals.h"
#define EPSZERO 1.e-8

/* read a double array a[n] into a Spamaline */ 
int darray2spamaline( struct Spamaline *sml, int n, double *a ) 
{
  int *jbuf, nsml, j, k ;

  jbuf = ( int * ) malloc( n * sizeof( int ) ) ;

  nsml = 0 ;
  for ( j = 0 ; j < n ; j++ )
    if ( fabs( a[j] ) > EPSZERO )
      jbuf[nsml++] = j ;
  sml -> n = nsml ;
  sml -> j = ( int * ) malloc( nsml * sizeof( int ) ) ;
  sml -> v = ( double * ) malloc( nsml * sizeof( double ) ) ;
  for ( k = 0 ; k < nsml ; k++ ) {
    j = jbuf[k] ;
    sml -> j[k] = j ;
    sml -> v[k] = a[j] ;
  }
  free( jbuf ) ;
  return nsml ;
}


int spamatranspose( struct Spamaline **psm, int nlines, int ncolumns ) 
{
  int *nextk, *ibuf, i, j, k, kt, counter = 0 ;
  double *vbuf ;
  struct Spamaline *sm, *smtransposed ;
  
  nextk = ( int * ) calloc( nlines, sizeof( int ) ) ;
  ibuf  = ( int * ) malloc( nlines * sizeof( int ) ) ;
  vbuf  = ( double * ) malloc( nlines * sizeof( double ) ) ;
  smtransposed = ( struct Spamaline * ) 
    malloc( ncolumns * sizeof( struct Spamaline ) ) ;
  sm = *psm ;

  for ( j = 0 ; j < ncolumns ; j++ ) {
    kt = 0 ;
    for ( i = 0 ; i < nlines ; i++ ) {
      k = nextk[i] ;
      if ( k < sm[i].n && j == sm[i].j[k] ) {  /* bingo: the element (i,j) is non zero */
	ibuf[kt] = i ;
	vbuf[kt] = sm[i].v[k] ;
	kt++ ;
	nextk[i]++ ;
      }
    }
    smtransposed[j].n = kt ;
    smtransposed[j].j = ( int * ) malloc( kt * sizeof( int ) ) ;
    smtransposed[j].v = ( double * ) malloc( kt * sizeof( double ) ) ;
    for ( k = 0 ; k < kt ; k++ ) {
      smtransposed[j].j[k] = ibuf[k] ;
      smtransposed[j].v[k] = vbuf[k] ;
    }
    counter += kt ;
  }

  /* deallocate  sm  == *psm  completely */  
  for ( i = 0 ; i < nlines ; i++ ) {
    free( sm[i].j ) ;
    free( sm[i].v ) ;
  }
  free( *psm ) ;

  *psm = smtransposed ;
  
  /* now, on exit, *psm is the transposed sparse matrix */

  return counter ;  /* total number of elements in  *psm  on exit */
}


void spamalinewrite( struct Spamaline sml ) 
{
  int k ;
  for ( k = 0 ; k < sml.n ; k++ ) 
    printf(" %2d %7.4lf,", sml.j[k], sml.v[k] ) ;
  printf("\n") ;
}
/* former version with pointer */
/*
void spamalinewrite( struct Spamaline *sml ) 
{
  int k ;
  for ( k = 0 ; k < sml -> n ; k++ ) 
    printf(" %2d %7.4lf,", sml -> j[k], sml -> v[k] ) ;
  printf("\n") ;
}
*/

void spamadelete( struct Spamaline *sparsemat, int nlines ) 
{
  int i ;
  for ( i = 0 ; i < nlines ; i++ ) {
    free( sparsemat[i].j ) ;
    free( sparsemat[i].v ) ;
  }
  free( sparsemat ) ;
}

/* 
   sparsematrealmatele   returns matrix element  <bra|A|ket>  
   where the operator  A  is described by psmline0,
   <bra| = \sum_i^{nbrabasis} bravec[i] <brabasis[i]|,  
   |ket> = \sum_j^{nketbasis} ketvec[j] |ketbasis[j]>   
   For consistency we must have that 
   psmline0 (as an array)  has dimension  nbrabasis  and that the 
   indices j = psmline0[i].j[k] do not exceed  nketbasis.
   <bra|A|ket> = \sum_i^{nbrabasis} \sum_j^{nketbasis} 
                 bravec[i] <brabasis[i]|A|ketbasis[j]> ketvec[j] 
        = \sum_i^{nbrabasis} bravec[i] *
          \sum_k^{psmline0[i].n} psmline0[i].v[k] * ketvec[psmline0[i].j[k]]

   WARNING : ONLY VALID WHEN COEFFICIENTS bravec[i] ARE REAL 
*/ 
/* 
  fast version. here we denote the pointer psmline rather then psmline0 
*/
double spamarealmatele( int nketbasis, double *ketvec, int nbrabasis,
     double *bravec, struct Spamaline *psmline ) 
{
  int i, k ;
  double sumk, sum ;

  sum = 0. ;  
  for ( i = 0 ; i < nbrabasis ; i++ ) {
    sumk = 0. ;
    for ( k = 0 ; k < psmline -> n ; k++ ) {
      sumk  +=  psmline -> v[k]  *  ketvec[ psmline -> j[k] ] ;
    }
    sum  +=  (*bravec++) * sumk ;
    psmline++ ;
  }
  return sum ;
}

/* more careful version */
/*
double spamarealmatele( int nketbasis, double *ketvec, int nbrabasis,
     double *bravec, struct Spamaline *psmline0 ) 
{
  int brai, k ;
  double sumk, sum ;
  struct Spamaline *psml ;

  sum = 0. ;
  for ( brai = 0 ; brai < nbrabasis ; brai++ ) {
    psml = psmline0 + brai ;
    sumk = 0. ;
    for ( k = 0 ; k < psml -> n ; k++ ) {
      sumk  +=  psml -> v[k]  *  ketvec[ psml -> j[k] ] ;
    }
    sum  +=  bravec[brai] * sumk ;
  }
  return sum ;
}
*/





/*  safer version : check for arrays of dimension 0 - but don't be so scared */
/*
void spamadelete( struct Spamaline *psparsematline0, int nlines ) 
{
  int i ;
  if ( nlines > 0 ) {
    for ( i = 0 ; i < nlines ; i++ ) 
      if ( psparsematline0[i].n > 0 ) {
        free( psparsematline0[i].j ) ;
        free( psparsematline0[i].v ) ;
      }
    free( psparsematline0 ) ;
  }
}
*/

/* deallocate a Spamaline  
void spamalinedealloc( struct Spamaline *sml ) 
{
    free( sml -> j ) ;
    free( sml -> v ) ;
}
*/


/* well functioning but no longer needed transpose of SQUARE sparse matrix */
/*
int spamasquaretranspose( struct Spamaline **psm, int nlines ) 
{
  int *nextk, *ibuf, i, j, k, kt, counter = 0 ;
  double *vbuf ;
  struct Spamaline *sm, *smtransposed ;
  
  nextk = ( int * ) calloc( nlines, sizeof( int ) ) ;
  ibuf  = ( int * ) malloc( nlines * sizeof( int ) ) ;
  vbuf  = ( double * ) malloc( nlines * sizeof( double ) ) ;
  smtransposed = ( struct Spamaline * ) 
    malloc( nlines * sizeof( struct Spamaline ) ) ;
  sm = *psm ;

  for ( j = 0 ; j < nlines ; j++ ) {
    kt = 0 ;
    for ( i = 0 ; i < nlines ; i++ ) {
      k = nextk[i] ;
      if ( k < sm[i].n && j == sm[i].j[k] ) {  
	ibuf[kt] = i ;
	vbuf[kt] = sm[i].v[k] ;
	kt++ ;
	nextk[i]++ ;
      }
    }
    smtransposed[j].n = kt ;
    smtransposed[j].j = ( int * ) malloc( kt * sizeof( int ) ) ;
    smtransposed[j].v = ( double * ) malloc( kt * sizeof( double ) ) ;
    for ( k = 0 ; k < kt ; k++ ) {
      smtransposed[j].j[k] = ibuf[k] ;
      smtransposed[j].v[k] = vbuf[k] ;
    }
    counter += kt ;
  }

  for ( i = 0 ; i < nlines ; i++ ) {
    free( sm[i].j ) ;
    free( sm[i].v ) ;
  }
  free( *psm ) ;

  *psm = smtransposed ;
  
  return counter ;  
}
*/
#undef EPSZERO
