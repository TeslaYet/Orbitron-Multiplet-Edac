#include <math.h>
#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
/* #define PI 3.141592653589793238512808959406186204433 */
#define FOURPI 12.56637061435917295405123583762474481773

double plgndr( int l, int m, double x ) 
{ /* computes the associated Legendre polynomial P_l^m(x). 0<=m<=l, -1.<x<1. */
  double fact, pll, pmm, pmmp1, somx2 ;
  int i, ll ;
  if ( l < 0 || m < 0 || m > l || fabs(x) > 1. ) {
    printf("bad args in plgdr()\n") ; exit(1) ; 
  }
  pmm = 1. ;
  if ( m > 0 ) {
    somx2 = sqrt( (1.-x)*(1.+x) ) ;
    fact = 1. ;
    for ( i = 1 ; i <= m ; i++ ) {
      pmm *= -fact*somx2 ;
      fact += 2. ;
    }
  }
  if ( l == m ) 
    return pmm ;
  else {
    pmmp1 = x * (2*m+1) * pmm ;
    if ( l == m + 1 ) 
      return pmmp1 ;
    else {
      for ( ll = m + 2 ; ll <= l ; ll++ ) {
	pll = ( x * (2*ll-1) * pmmp1 - (ll+m-1) * pmm ) / ( ll - m ) ;
        pmm = pmmp1 ;
        pmmp1 = pll ;
      }
      return pll ;
    }
  }
}    
double complex cylm( int l, int m, double x, double phi)
{
  int signm, absm, k, sign ;
  double dum, yre, yim ;
  if ( l < 0 || m < -l || m > l || fabs(x) > 1. ) {
    printf("bad args in cylm()\n") ; exit(1) ; 
  }
  signm = ( m >= 0 ? 1 : -1 ) ;
  absm  = signm * m ;
  dum = 1. ;
  for ( k = l - absm + 1 ; k <= l + absm ; k++ )
    dum *= k ;  
  dum = sqrt( ( 2. * l + 1. ) / FOURPI / dum ) ;
  dum *= plgndr( l, absm, x ) ;
  yre = dum * cos( absm * phi ) ;
  yim = dum * sin( absm * phi ) ;
  if ( signm == -1 ) {
    sign = ( absm % 2 == 0 ? 1 : -1 ) ;
    yre *=  sign ;
    yim *= -sign ;   
  }
  return yre + I*yim ;
}

