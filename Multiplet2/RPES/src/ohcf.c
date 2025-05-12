/* applies crystal field in Oh symmetry */
#include <stdlib.h>
#include <stdio.h>
#include "globals.h"
int crystalfieldOh2list( int nshells, int *lsh, int *sorb1sh, 
			 double *tendqsh, struct O1plistitem **ppo1plistitem0 ) 
{
  int ish, l0, m, sigma, i1, i2, counter = 0 ;
  double tendq, shift ;

  for ( ish = 0 ; ish < nshells ; ish++ ) {
    l0   = lsh[ish] ; 
    tendq = tendqsh[ish] ;
    if ( tendq != 0. && l0 > 1 ) {
      if ( l0 == 2 ) {
      /* diagonal elements: 0:eg. 1,-1:t2g,  <2|CF|2>=<-2|CF|-2>=(eg+t2g)/2 */
        for ( m = -l0 ; m <= l0 ; m++ ) {
	  if ( m == 0 ) shift = 0.6 * tendq ;
	  else if ( m == 1 || m == -1 ) shift = -0.4 * tendq ;
	  else if ( m == 2 || m == -2 ) shift =  0.1 * tendq ;
	  else { printf("TROUBLE in crystalfieldOh2list\n") ; exit(1) ; }
	  for ( sigma = -1 ; sigma <= 1 ; sigma += 2 ) {
	    i1 = sporbindex( sorb1sh, lsh, ish, m, sigma ) ; 
	    o1plistitemadd( ppo1plistitem0, i1, i1, shift ) ;
	    counter++ ;
	  }      
	}
	/* off diagonal elements: <2|CF|-2> = <-2|CF|2> = (eg-t2g)/2 */
	shift = 0.5 * tendq ;
	for ( sigma = -1 ; sigma <= 1 ; sigma += 2 ) {
	  i1 = sporbindex( sorb1sh, lsh, ish, 2, sigma ) ;
	  i2 = sporbindex( sorb1sh, lsh, ish,-2, sigma ) ;
	  o1plistitemadd( ppo1plistitem0, i1, i2, shift ) ;
	  o1plistitemadd( ppo1plistitem0, i2, i1, shift ) ;
	  counter += 2 ;
	}      
      }
      else { /* l > 2 */
	printf("CF not applied to shell No %d. Don't know how for l = %d\n",
	       ish, l0  ) ;
      }
    }
  }
  return counter ;
}
