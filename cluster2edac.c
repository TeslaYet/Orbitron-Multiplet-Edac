#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define SYMBOLLENGTH 5
#define FACTOR 5
#define EPSILON 1.e-6
#define BOHR 0.529177

struct Atom { char t[SYMBOLLENGTH] ; double x, y, z ; } ;

double veclength( double vec[3] ) ;

double distance( struct Atom atom1, struct Atom atom2) 
{
    double dvec[3] ;
    dvec[0] = atom1.x - atom2.x ;
    dvec[1] = atom1.y - atom2.y ;
    dvec[2] = atom1.z - atom2.z ;
    return veclength( dvec ) ;
}

void calcastar( double a[3][3], double as[3][3] ) ;

void calcsurfnormal( double as[3][3], int miller[3], double surfnormal[3]);

void rotate( double cosph, double sinph, double costh, double sinth,
	     double vec[3] ) ;

void sortatoms( int natoms, struct Atom *atom ) ;

int main()
{
  char **type ;
  int nbasisatoms, i, j, k, miller[3], natoms = 0, n = 0 ;
  int i0,i1,i2,i0max,i1max,i2max ;
  double a[3][3], as[3][3], ap[3][3], lengtha[3], **red ;
  double sum, surfnormal[3], y[3], cutoff, radius2, depth2, r2, zsurf ;
  double originred[3], origin[3] ;
  double  latticeconstant, dum, cosph, sinph, costh, sinth ;
  struct Atom *atom ;
  FILE *fp ;

  for ( i=0 ; i<3 ; i++ ) {
    printf("%dth lattice vector a%d_x, a%d_y, a%d_z: ",i,i,i,i) ;
    for ( j=0 ; j<3 ; j++ )
      scanf("%lf",&a[i][j]) ;
  }
  calcastar(a,as);
  for ( i=0 ; i<3 ; i++ ) {
    printf("a_%d = ", i ) ;
    for ( j=0 ; j<3 ; j++ )
      printf("%7.3lf ",a[i][j]) ;
    printf(" \ta*_%d = ", i ) ;
    for ( j=0 ; j<3 ; j++ ) 
      printf("%7.3lf ",as[i][j]) ;
    printf("\n");
  }

  printf("\n");
  printf("Miller indices of surface (hkl): ") ; 
  for ( i = 0 ; i < 3 ; i++ ) 
    scanf("%d", &miller[i] ) ;

  calcsurfnormal( as, miller, surfnormal ) ;

  /* rotation such that surfnormal -> e_z */
  dum = sqrt( surfnormal[0]*surfnormal[0] + surfnormal[1]*surfnormal[1] );
  if ( dum > 1.e-10 ) {
    cosph = surfnormal[0] / dum ; 
    sinph = surfnormal[1] / dum ; 
    costh = surfnormal[2] / veclength( surfnormal ) ;
    sinth = dum / veclength( surfnormal ) ;
  } else {
    costh = 1. ;
    sinth = 0. ;
    cosph = 1. ;
    sinph = 0. ;
  }
  /*  printf("1-cosph^2+sinph^2 = %lf\n", 1.-cosph*cosph-sinph*sinph );
  printf("1-costh^2+sinth^2 = %lf\n", 1.-costh*costh-sinth*sinth );
  */
  for ( i = 0 ; i < 3 ; i++ )
    for ( j = 0 ; j < 3 ; j++ )
      ap[i][j] = a[i][j] ;
  for ( i = 0 ; i < 3 ; i++ ) 
    rotate( cosph, sinph, costh, sinth, ap[i] ) ;

  for ( i=0 ; i<3 ; i++ ) {
    printf("a_%d = ", i ) ;
    for ( j=0 ; j<3 ; j++ )
      printf("%7.3lf ",a[i][j]) ;
    printf(" \ta'_%d = ", i ) ;
    for ( j=0 ; j<3 ; j++ ) 
      printf("%7.3lf ",ap[i][j]) ;
    printf("\n");
  }

  printf("\norigin of frame in red. coords? " ) ;
  for ( i = 0 ; i < 3 ; i++ ) 
    scanf("%lf", &originred[i] ) ;
  printf("\norigin in cart. coords = ") ;
  for ( j = 0 ; j < 3; j++ ) {
    sum = 0. ;
    for ( i = 0 ; i < 3 ; i++ )
      sum += a[i][j] * originred[i] ;
    printf(" %7.3lf", sum ) ;
  }
  printf("\n") ;

  printf("\nNo of basis atoms: ") ; scanf("%d",&nbasisatoms) ;
  type = ( char ** ) malloc( nbasisatoms * sizeof( char * ) ) ;
  red  = ( double ** ) malloc( nbasisatoms * sizeof( double * ) ) ;
  for ( i=0 ; i<nbasisatoms ; i++ ) {
    type[i] = ( char * ) malloc( (SYMBOLLENGTH+1) * sizeof( char ) ) ;
    red[i] = ( double * ) malloc( 3 * sizeof( double ) ) ;
    printf("\nSymbol and reduced coords of atom No %d: ", i ) ;
    scanf("%s%lf%lf%lf", type[i], &red[i][0], &red[i][1], &red[i][2] ) ;
    printf("%s %lf %lf %lf\n", type[i], red[i][0], red[i][1], red[i][2] ) ;
  }
  printf("Positions in cartesian coordinates\n") ;
  for ( k=0 ; k<nbasisatoms ; k++ ) {
    printf("\n%s", type[k] ) ;
    for ( j = 0 ; j < 3 ; j++ ) {
      sum = 0. ;
      for ( i = 0 ; i < 3 ; i++ )
        sum += a[i][j] * red[k][i] ;
      printf("\t%lf", sum ) ;
    }
  }

  printf("lattice constant in Angs? ") ; scanf("%lf",&latticeconstant) ;
  printf("\nradius of ellipsoid in Angs? ") ; scanf("%lf", &radius2 ) ;
  cutoff = radius2 ;
  radius2 *= radius2 ;
  printf("\ndepth of ellipsoid in Angs? ") ; scanf("%lf", &depth2 ) ;
  if ( depth2 > cutoff ) cutoff = depth2 ;
  depth2 *= depth2 ;
  printf("\nzsurf in Angs? ") ; scanf("%lf", &zsurf ) ;

  i0max = cutoff * FACTOR / veclength( a[0] ) ;
  i1max = cutoff * FACTOR / veclength( a[1] ) ;
  i2max = cutoff * FACTOR / veclength( a[2] ) ;

  printf("\nCoordinates\n") ;
  for ( i0 = -i0max ; i0 <= i0max ; i0++ ) {
    for ( i1 = -i1max ; i1 <= i1max ; i1++ ) {
      for ( i2 = -i2max ; i2 <= i2max ; i2++ ) {
        for ( k = 0 ; k < nbasisatoms ; k++ ) {
          for ( j = 0 ; j < 3 ; j++ ) {
            y[j] =  ( i0 + red[k][0] - originred[0] ) * ap[0][j] 
                  + ( i1 + red[k][1] - originred[1] ) * ap[1][j]
                  + ( i2 + red[k][2] - originred[2] ) * ap[2][j] ;
            y[j] *= latticeconstant ;
          }
          r2 = ( y[0]*y[0] + y[1]*y[1] ) / radius2 + y[2]*y[2] / depth2 ;
          if ( r2 < 1. && y[2] < zsurf ) {
            printf("%4s %6.3lf %6.3lf %6.3lf\n", type[k],y[0],y[1],y[2] ) ;
	    natoms++ ;
	  }
        }
      }
    }
  }
  printf("\n");

  
  atom = ( struct Atom * ) malloc( natoms * sizeof( struct Atom ) ) ;
  n = 0 ;
  for ( i0 = -i0max ; i0 <= i0max ; i0++ ) {
    for ( i1 = -i1max ; i1 <= i1max ; i1++ ) {
      for ( i2 = -i2max ; i2 <= i2max ; i2++ ) {
        for ( k = 0 ; k < nbasisatoms ; k++ ) {
          for ( j = 0 ; j < 3 ; j++ ) {
            y[j] =  ( i0 + red[k][0] - originred[0] ) * ap[0][j] 
                  + ( i1 + red[k][1] - originred[1] ) * ap[1][j]
                  + ( i2 + red[k][2] - originred[2] ) * ap[2][j] ;
            y[j] *= latticeconstant ;
          }
          r2 = ( y[0]*y[0] + y[1]*y[1] ) / radius2 + y[2]*y[2] / depth2 ;
          if ( r2 < 1. && y[2] < zsurf ) {
            strcpy( atom[n].t, type[k] ) ; 
	    atom[n].x = y[0] ;
	    atom[n].y = y[1] ;
	    atom[n].z = y[2] ;
	    n++ ;
	  }
        }
      }
    }
  }
  /* printf("natoms = %d, n = %d\n", natoms, n ) ; */

  fp = fopen("cluster.xyz","w") ;
  fprintf(fp," %d\n\n", natoms);
  sortatoms( natoms, atom ) ;
  for ( i = 0 ; i < natoms ; i++ )
      fprintf(fp,"%4s %10.6lf %10.6lf %10.6lf\n",
	   atom[i].t, atom[i].x, atom[i].y, atom[i].z ) ; 
  fclose(fp);
  fp = fopen("xxx.clus","w") ;
  fprintf(fp,"%4d l(A)\n", natoms ) ;
  for ( i = 0 ; i < natoms ; i++ )
      fprintf(fp,"%3d %4s %13.9lf %13.9lf %13.9lf\n",
	      i+1, atom[i].t, atom[i].x, atom[i].y, atom[i].z ) ;
  fclose(fp);
  return 0 ;
}

void calcastar( double a[3][3], double as[3][3] ) 
{
  int i, j ;
  double volume ;
  for ( i = 0 ; i < 3 ; i++ ) 
    for ( j = 0 ; j < 3 ; j++ ) 
      as[i][j] = a[(i+1)%3][(j+1)%3] * a[(i+2)%3][(j+2)%3]
	- a[(i+1)%3][(j+2)%3] * a[(i+2)%3][(j+1)%3] ;
  volume = 0. ;
  for ( j = 0 ; j < 3 ; j++ )
    volume += a[0][j]*as[0][j] ;
  printf("cell volume = %lf\n", volume ) ;
  for ( i = 0 ; i < 3 ; i++ )
    for ( j = 0 ; j < 3 ; j++ )
      as[i][j] /= volume ;
}

void calcsurfnormal( double as[3][3], int miller[3], double surfnormal[3] )
{
  int i, j ;
  double dum ;
  for ( j = 0 ; j < 3; j++ ) {
    dum = 0. ;
    for ( i = 0 ; i < 3 ; i++ )
      dum += as[i][j] * miller[i] ;
    surfnormal[j] = dum ;
  }
  dum = veclength( surfnormal ) ;
  for ( j = 0 ; j < 3; j++ ) 
    surfnormal[j] /= dum ;
  printf("\nsurface normal = ") ;
  for ( j = 0 ; j < 3; j++ ) printf(" %7.3lf", surfnormal[j]) ;
  printf("\n") ;
}

void rotate( double cosph, double sinph, double costh, double sinth,
	     double vec[3] ) 
{
  int j ;
  double dum[3] ;
  for ( j = 0 ; j < 3 ; j++ ) dum[j] = vec[j] ;
  vec[0] = costh*cosph * dum[0] + costh*sinph * dum[1] - sinth * dum[2] ;
  vec[1] =      -sinph * dum[0] +       cosph * dum[1] ; 
  vec[2] = sinth*cosph * dum[0] + sinth*sinph * dum[1] + costh * dum[2] ;
}

double veclength( double vec[3] ) {
  int j ;
  double sum = 0. ;
  for ( j = 0 ; j < 3; j++ ) 
      sum += vec[j]*vec[j] ;
  return sqrt( sum ) ;
} 

/* bubble sort */
int atom1before2( struct Atom a, struct Atom b ) {
  return ( a.x*a.x+a.y*a.y+1.001*a.z*a.z < b.x*b.x+b.y*b.y+1.001*b.z*b.z ) ;
}
/*
int atom1before2( struct Atom a, struct Atom b ) {
  return ( a.x*a.x+a.y*a.y+a.z*a.z < b.x*b.x+b.y*b.y+b.z*b.z ) ;
} 
*/
/*
int atom1before2( struct Atom a, struct Atom b ) {
  return ( a.z > b.z + EPSILON || ( a.z > b.z - EPSILON && 
	   a.x*a.x+a.y*a.y < b.x*b.x+b.y*b.y ) ) ;
}
*/
 
void sortatoms( int natoms, struct Atom *atom )
{
  int sorted = 0, i ;
  struct Atom dumatom ;
  while ( !sorted ) {
    sorted = 1 ;
    for ( i = 0 ; i < natoms-1 ; i++ ) {
      if ( atom1before2( atom[i+1], atom[i] ) ) {
	sorted = 0 ;
	dumatom = atom[i] ;
	atom[i] = atom[i+1] ;
	atom[i+1] = dumatom ;
      }
    }
  }
}

/*
void rotate( double cosph, double sinph, double costh, double sinth,
	     double vec[3] ) 
{
  int j ;
  double duma[3], dumb[3] ;
  for ( j = 0 ; j < 3 ; j++ ) duma[j] = vec[j] ;
  dumb[0] =  cosph * duma[0] + sinph * duma[1] ;
  dumb[1] = -sinph * duma[0] + cosph * duma[1] ;
  dumb[2] =  duma[2] ;
  vec[0]  =  costh * dumb[0] - sinth * dumb[2] ;
  vec[1]  =  dumb[1] ;
  vec[2]  =  sinth * dumb[0] + costh * dumb[2] ;
}
*/
/*
void rotzminusph( double cosph, double sinph, double vec[3] )
{
  int j ;
  double dum[3] ;
  for ( j = 0 ; j < 3 ; j++ ) 
    dum[j] = vec[j] ;
  vec[0] =  cosph * dum[0] + sinph * dum[1] ;
  vec[1] = -sinph * dum[0] + cosph * dum[1] ;
  vec[2] =  dum[2] ;
}
 
void rotyminusth( double costh, double sinth, double vec[3] )
{
  int j ;
  double dum[3] ;
  for ( j = 0 ; j < 3 ; j++ ) 
    dum[j] = vec[j] ;
  vec[0] =  costh * dum[0] - sinth * dum[2] ;
  vec[1] =  dum[1] ;
  vec[2] =  sinth * dum[0] + costh * dum[2] ;
}
*/
