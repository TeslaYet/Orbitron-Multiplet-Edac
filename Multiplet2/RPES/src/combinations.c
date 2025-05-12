#include <stdlib.h>
/*********************************************************************/
/*   noverk returns the combinatorial factor ( n over k )            */ 
/*   = n! / ( k! (n-k)! )   and   0   if  n < 0 || k < 0 || k > n    */
/*********************************************************************/
int noverk( int n, int k )
{
  int res, j, m ;
  if ( n < 0 || k < 0 || k > n ) 
    return 0 ;
 
  res = 1 ; j = 1 ; m = n ;
  while ( j <= k ) {
    res *= m-- ;
    res /= j++ ;
  }
  return res ;
}

/*********************************************************************/
/*   combinations calculates all (n over k) combinations to choose   */
/*   k numbers from a total of n. The n numbers to choose from are   */
/*   offset, offset+1, ... , offset + n - 1.                         */
/*   result stored in 2-dim array  comb[i][j], i=0,noverk-1; j=0,k-1 */
/*********************************************************************/
int combinations( int n, int k, int offset, int ***comb ) 
{  
  int nok, i, j, jrun, *l ;
   
  nok  = noverk( n, k ) ;

  if ( nok > 0 ) {
    /* memory allocation */
    *comb = ( int ** ) malloc( nok * sizeof( int * ) ) ;
    for ( i = 0 ; i < nok ; i++ )
      (*comb)[i] = ( int * ) malloc( k * sizeof( int ) ) ; 

    l = ( int * ) malloc( k * sizeof( int ) ) ;
    /* initialization */
    for ( j = 0 ; j < k ; j++ )
      l[j] = j ;
    jrun = k - 1 ;
    /* main loop */
    i = 0 ;
    while ( jrun >= 0 ) {
      for ( j = 0 ; j < k ; j++ ) 
	(*comb)[i][j] = l[j] + offset ;
      while ( jrun >= 0 && l[jrun] == n - k + jrun ) 
	jrun-- ;
      if ( jrun >= 0 ) {
	l[jrun]++ ;
	for ( j = jrun + 1 ; j < k ; j++ )
	  l[j] = l[j-1] + 1 ;
	jrun = k - 1 ;
      }
      i++ ;
    }
    /* memory deallocation */
    free( l ) ;
  }
  return nok ;
}
