#include <stdio.h>
#include <math.h>
#include "w3j.h"
#include "globals.h"

/* calculate \int_0^r_0 P_inish(r) r P_finsh(r) dr */
/*double radipmatele( int inish, int finsh ) 
{
  return 1. ;
}
*/
/* 
   dipmatele2list  calculates the matrix elements of the dipole operator 
   r^(1)_q == r C^(1)_q, using formulae (11.15), (14.58) of Cowan app F-7. 
   for given  nli  and q  (where  l,m <=> lf,mf and l',m' <=> li,mi).
   It writes the non-zero elements to a linked list. 
*/

int dipole2list( int nshells, int *lsh, int *sorb1sh, int inish, int q,
                 double **radipmatele, struct O1plistitem **ppitem0 )
{
  int finsh, li, lf, lgreater, sign, mi, mf, sigma, iorbi, iorbf ;
  int counter = 0 ;
  double reducedmatele, val ;
  
  /* uses formulae (14.58) (11.15) of Cowan app F.7 */
  /* with l,m <=> lf,mf and l',m' <=> li,mi         */

  li = lsh[inish] ;
  for ( finsh = 0 ; finsh < nshells ; finsh++ ) {
    lf = lsh[finsh] ;
    if ( lf == li + 1 || lf == li - 1 ) {
      /* calculate reduced matrix element using (14.58) */
      lgreater = li > lf ? li : lf ;
      sign =  ( lf + lgreater ) % 2 == 0 ? 1 : -1 ;
      reducedmatele = sign * sqrt( lgreater ) * radipmatele[inish][finsh] ;
      /* calculate matrix element using Wigner Eckart Theorem (11.15) */
      for (  mi = -li ; mi <= li ; mi++ ) {
	mf = mi + q ;
	if ( -lf <= mf && mf <= lf ) {
	  sign = ( lf - mf ) % 2 == 0 ? 1 : -1 ;
	  val = sign * w3j( 2*lf, 2, 2*li, -2*mf, 2*q, 2*mi ) * reducedmatele ;
	printf("li=%d, mi=%d, lf=%d, mf=%d, si=%d, w3j=%lf, re=%lf, val=%lf\n",
	li,mi,lf,mf,sign,w3j(2*lf,2,2*li,-2*mf,2*q,2*mi),reducedmatele,val ) ;
	  if ( val != 0 ) {
	    /* write matrix elements in list */
	    for ( sigma = -1 ; sigma <= 1 ; sigma += 2 ) {
	      iorbi = sporbindex( sorb1sh, lsh, inish, mi, sigma ) ;
	      iorbf = sporbindex( sorb1sh, lsh, finsh, mf, sigma ) ;
	      o1plistitemadd( ppitem0, iorbi, iorbf, val ) ;
	      counter ++ ;
	    }
	  }
	}
      }
    }
  }
  return counter ;
}
