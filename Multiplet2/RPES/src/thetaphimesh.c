#include <stdlib.h>
#include <math.h>
#define RADEG .01745329251994329576
#define EPSILON 1.e-5
void thetaphimesh( double thmin, double thmax, double thdelta, 
                   double phmin, double phmax, double phdelta, 
                int *pntheta, double **ptheta, int **pnphi, double ***pphi ) 
{
  int ntheta, *nphi, i, j ;
  double phde, *theta, **phi ; 
  thmax *= RADEG ; thmin *= RADEG ; thdelta *= RADEG ;
  phmax *= RADEG ; phmin *= RADEG ; phdelta *= RADEG ;
  ntheta = ( int ) ( (thmax-thmin)/thdelta + 1.5 ) ;
  if ( ntheta > 1 ) thdelta = ( thmax - thmin) / ( ntheta - 1 ) ; 
  theta = ( double * ) malloc( ntheta*sizeof( double ) ) ;
  nphi  = ( int * ) malloc( ntheta*sizeof( int ) ) ;
  phi = ( double ** ) malloc( ntheta*sizeof( double * ) ) ;
  for ( i = 0 ; i < ntheta ; i++ ) {
    theta[i] = thmin + i * thdelta ;
    if ( theta[i] < EPSILON ) { phde = 0. ; nphi[i] = 1 ; } 
    else { 
      phde = phdelta/sin(theta[i]) ;
      nphi[i] = ( int ) ( (phmax-phmin)/phde + 1.5 )  ;
      if ( nphi[i] > 1 ) phde = ( phmax - phmin) / ( nphi[i] - 1 ) ;
    }
    phi[i] = ( double * ) malloc( nphi[i]*sizeof( double ) ) ;
    for ( j = 0 ; j < nphi[i] ; j++ )  phi[i][j] = phmin + j * phde ;
  }
  *pntheta = ntheta ;
  *ptheta = theta ;
  *pnphi = nphi ;
  *pphi = phi ; 
}
#undef EPSILON
