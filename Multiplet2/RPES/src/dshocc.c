/* counts d-orbitals of some symmetry */
#include <stdlib.h>
#include <stdio.h>
#include "globals.h"
int dxyocc2list( int nshells, int *lsh, int *sorb1sh, int ish,
                         struct O1plistitem **ppo1plistitem0 ) 
{
  int m1, m2, sigma, i1, i2, counter = 0 ;
  double coeff[5][5] = { 
        { 0.5, 0., 0., 0.,-0.5 },
        { 0. , 0., 0., 0., 0.  },
        { 0. , 0., 0., 0., 0.  },
        { 0. , 0., 0., 0., 0.  },
        {-0.5, 0., 0., 0., 0.5 } } ;

  if ( lsh[ish] != 2 ) { 
     printf("dxyocc2list: Shell No %d is not a d-shell. Stop.\n",ish);
     exit(1) ;
  }
  for ( m1 = -2 ; m1 <= 2 ; m1++ ) {
    for ( m2 = -2 ; m2 <= 2 ; m2++ ) {
      for ( sigma = -1 ; sigma <= 1 ; sigma += 2 ) {
        i1 = sporbindex( sorb1sh, lsh, ish, m1, sigma ) ;
        i2 = sporbindex( sorb1sh, lsh, ish, m2, sigma ) ;
        if ( coeff[m1+2][m2+2] != 0. ) {
          o1plistitemadd( ppo1plistitem0, i1, i2, coeff[m1+2][m2+2] ) ; 
          counter++ ;
	}
      }
    }
  }
  return counter ;
}
int t2gocc2list( int nshells, int *lsh, int *sorb1sh, int ish,
                         struct O1plistitem **ppo1plistitem0 ) 
{
  int m1, m2, sigma, i1, i2, counter = 0 ;
  double coeff[5][5] = { 
        { 0.5, 0., 0., 0.,-0.5 },
        { 0. , 1., 0., 0., 0.  },
        { 0. , 0., 0., 0., 0.  },
        { 0. , 0., 0., 1., 0.  },
        {-0.5, 0., 0., 0., 0.5 } } ;

  if ( lsh[ish] != 2 ) { 
     printf("t2gocc2list: Shell No %d is not a d-shell. Stop.\n",ish);
     exit(1) ;
  }
  for ( m1 = -2 ; m1 <= 2 ; m1++ ) {
    for ( m2 = -2 ; m2 <= 2 ; m2++ ) {
      for ( sigma = -1 ; sigma <= 1 ; sigma += 2 ) {
        i1 = sporbindex( sorb1sh, lsh, ish, m1, sigma ) ;
        i2 = sporbindex( sorb1sh, lsh, ish, m2, sigma ) ;
        if ( coeff[m1+2][m2+2] != 0. ) {
          o1plistitemadd( ppo1plistitem0, i1, i2, coeff[m1+2][m2+2] ) ; 
          counter++ ;
	}
      }
    }
  }
  return counter ;
}
int egocc2list( int nshells, int *lsh, int *sorb1sh, int ish,
                         struct O1plistitem **ppo1plistitem0 ) 
{
  int m1, m2, sigma, i1, i2, counter = 0 ;
  double coeff[5][5] = { 
        { 0.5, 0., 0., 0., 0.5 },
        { 0. , 0., 0., 0., 0.  },
        { 0. , 0., 1., 0., 0.  },
        { 0. , 0., 0., 0., 0.  },
        { 0.5, 0., 0., 0., 0.5 } } ;

  if ( lsh[ish] != 2 ) { 
     printf("t2gocc2list: Shell No %d is not a d-shell. Stop.\n",ish);
     exit(1) ;
  }
  for ( m1 = -2 ; m1 <= 2 ; m1++ ) {
    for ( m2 = -2 ; m2 <= 2 ; m2++ ) {
      for ( sigma = -1 ; sigma <= 1 ; sigma += 2 ) {
        i1 = sporbindex( sorb1sh, lsh, ish, m1, sigma ) ;
        i2 = sporbindex( sorb1sh, lsh, ish, m2, sigma ) ;
        if ( coeff[m1+2][m2+2] != 0. ) {
          o1plistitemadd( ppo1plistitem0, i1, i2, coeff[m1+2][m2+2] ) ; 
          counter++ ;
	}
      }
    }
  }
  return counter ;
}
