/* applies crystal field to d-shell */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "globals.h"
#define EPSILON 1.e-7
/* CF as real 5x5 matrix cfdmat. Trace is subtracted */
int crystalfield2list( int nshells, int *lsh, int *sorb1sh, int ish,
		       double **cfdmat, struct O1plistitem **ppo1plistitem0 ) 
{
  int m1, m2, sigma, i1, i2, counter = 0 ;
  double dum ;

  if ( lsh[ish] != 2 ) { 
    printf("crystalfield2list: Shell No %d is not a d-shell. Stop.\n",ish);
    exit(1) ;
  }
  /* make cfdmat traceless */
  for ( dum = 0., m1 = -2 ; m1 <= 2 ; m1++)
    dum += cfdmat[2+m1][2+m1] ;
  for ( m1 = -2 ; m1 <= 2 ; m1++)
    cfdmat[2+m1][2+m1] -= dum / 5. ;

  for ( m1 = -2 ; m1 <= 2 ; m1++) {     
    for ( m2 = -2 ; m2 <= 2 ; m2++) {
      dum = cfdmat[2+m1][2+m2] ;
      if ( fabs( dum ) > EPSILON ) 
	for ( sigma = -1 ; sigma <= 1 ; sigma += 2 ) {
	  i1 = sporbindex( sorb1sh, lsh, ish, m1, sigma ) ; 
	  i2 = sporbindex( sorb1sh, lsh, ish, m2, sigma ) ; 
	  o1plistitemadd( ppo1plistitem0, i1, i2, dum ) ; 
	  counter++ ;
	}
    }
  }
  return counter ;
}
/* d-shell CF for C4 symmetry or higher */
int crystalfieldC42list( int nshells, int *lsh, int *sorb1sh, int ish,
                         double ea1, double eb1, double eb2, double ee,
                         struct O1plistitem **ppo1plistitem0 ) 
{
  int m, sigma, i1, i2, counter = 0 ;
  double shift ;

  if ( lsh[ish] != 2 ) { 
     printf("crystalfieldC42list: Shell No %d is not a d-shell. Stop.\n",ish);
     exit(1) ;
  }
  shift = 0.2*(ea1+eb1+eb2+2.*ee) ;
  ea1 -= shift ; eb1 -= shift ; eb2 -= shift ; ee -= shift ;
  printf("Eav substracted. a1=%.4lf b1=%.4lf b2=%.4lf e=%.4lf\n", ea1,eb1,eb2,ee) ;
  /* diagonal elements: 0:ea1, 1,-1:ee,  <2|CF|2>=<-2|CF|-2>=(eb1+eb2)/2 */
  for ( m = -2 ; m <= 2 ; m++ ) {
    if ( m == 0 ) shift = ea1 ;
    else if ( m == 1 || m == -1 ) shift = ee ;
    else if ( m == 2 || m == -2 ) shift = 0.5*(eb1+eb2) ;
    else { printf("TROUBLE in crystalfieldC42list\n") ; exit(1) ; }
    if ( fabs( shift ) > EPSILON ) {
      for ( sigma = -1 ; sigma <= 1 ; sigma += 2 ) {
        i1 = sporbindex( sorb1sh, lsh, ish, m, sigma ) ; 
        o1plistitemadd( ppo1plistitem0, i1, i1, shift ) ;
        counter++ ;
      }      
    }
  }
  /* off diagonal elements: <2|CF|-2> = <-2|CF|2> = (eb1-eb2)/2 */
  shift = 0.5*(eb1-eb2) ;
  if ( fabs( shift ) > EPSILON ) { 
    for ( sigma = -1 ; sigma <= 1 ; sigma += 2 ) {
      i1 = sporbindex( sorb1sh, lsh, ish, 2, sigma ) ;
      i2 = sporbindex( sorb1sh, lsh, ish,-2, sigma ) ;
      o1plistitemadd( ppo1plistitem0, i1, i2, shift ) ;
      o1plistitemadd( ppo1plistitem0, i2, i1, shift ) ;
      counter += 2 ;
    }      
  }
  return counter ;
}
#undef EPSILON
