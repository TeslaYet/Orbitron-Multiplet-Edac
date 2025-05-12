#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define MARGIN 0
#define THETAFRAME 82
#define POINTSPERDEG 4
#define EPSILON 1.e-8
#define PI 3.141592653589793116
#define LINECHMAX 120

double interpol( double theta, double phi, int ntheta, double dtheta,
                 int nphi, double dphi, double **chi ) ;

void writechiline(FILE *fp, int chilinelength, int *chiline ) ;

int main( int argc, char *argv[] )
{
  char pngfile[30], command[80] ;
  char col16[16][3]={ "00","11","22","33","44","55","66","77","88","99",
		      "aa","bb","cc","dd","ee","ff" } ;
  int ntheta, nphi, i, j, counter, itheta, jphi, np ;
  int npoints, chilinelength, *chiline, ipt, jpt ;
  int iy80, imy80, jx80, jmx80, ichichi ;
  double thetamax, dtheta, dphi, theta, phi, intensin, phimax, theta0, phi0 ;
  double  **chi, chichi, th, ph, x, y, r ;
  double chiij, chimin, chimax, chiav, chivar ;
  double chiminin, chimaxin ;
  char line[LINECHMAX+1] ;
  FILE *fp ;

  fp = fopen( (argc > 1?argv[1]:"output.ms"),"r") ; 
  fgets(line,LINECHMAX,fp) ; fgets(line,LINECHMAX,fp) ;
  thetamax = phimax = 0.;
  dtheta = dphi = 0.;
  np = 0 ;
  fscanf(fp,"%lf%lf%lf",&theta0,&phi0,&intensin) ;
  np++ ;
  fscanf(fp,"%lf%lf%lf",&theta,&phi,&intensin) ;
  np++ ;
  if ( theta > theta0 ) dtheta = theta - theta0 ;
  else if ( phi > phi0 ) dphi = phi - phi0 ;
  else { printf("Neither theta nor phi advance. STOP\n") ; exit(1) ; }
  
  while( fscanf(fp,"%lf%lf%lf",&theta,&phi,&intensin) != EOF ) {
    if ( theta > thetamax ) thetamax = theta ;
    if ( phi > phimax ) phimax = phi ;
    np++ ; 
  }
  fclose(fp) ;

  if ( dphi > 0. ) {
    nphi = (int) (phimax-phi0)/dphi + 1.5 ;
    ntheta = np / nphi ;
    dtheta = (thetamax-theta0)/(ntheta-1) ;
  } else {
    ntheta = (int) (thetamax-theta0)/dtheta + 1.5 ;
    nphi = np / ntheta ;
    dphi = (phimax-phi0)/(nphi-1) ;
  }
  if ( fabs(phimax-360.) > EPSILON ||
       fabs(theta0) > EPSILON ||
       fabs(phi0) > EPSILON ) {
    printf("can't use this plot program. STOP\n") ;
    exit(1) ;
  }
  /*  printf("%lf %lf %lf %lf %d\n",dtheta,dphi,thetamax,phimax,np) ;
  printf("%d %d\n",ntheta,nphi);
  exit(1) ; */
  
  chiminin = 1. ; chimaxin= -1. ; /* automatic min/max */

  printf("%lf %lf %lf %d %d\n", thetamax, dtheta, dphi, ntheta, nphi ) ;

  chi = ( double ** ) calloc( ntheta, sizeof( double * ) ) ;
  for ( i = 0 ; i < ntheta ; i++ )
    chi[i] = ( double * ) calloc( nphi, sizeof( double ) ) ;

  fp = fopen( (argc > 1?argv[1]:"output.ms"),"r") ; 
  fgets(line,LINECHMAX,fp) ; fgets(line,LINECHMAX,fp) ;
  for ( i = 0 ; i < ntheta ; i++ ) {
    for ( j = 0 ; j < nphi ; j++ ) {
      fscanf(fp,"%lf%lf%lf", &theta, &phi, &chi[i][j] ) ;
      if ( fabs(i * dtheta - theta ) > EPSILON ||
           fabs(j * dphi - phi ) > EPSILON ) {printf("i=%d j=%d.CHAOS\n",i,j);exit(1);}
    }
    /* chi[i][nphi-1] = chi[i][0] ; DON't put phi==360 as phi==0 */
  }
  fclose(fp) ;
  fp = fopen("ordered.dat","w") ;
  fprintf(fp,"%4lg %4lg  %lg\n", thetamax, dtheta, dphi ) ;
  for ( i = 1 ; i < ntheta ; i++ ) {
    itheta = i*dtheta + EPSILON ;
    for ( j = 0 ; j < nphi - 1 ; j++ ) {
      jphi = j*dphi + EPSILON ;
      fprintf(fp,"%4d %4d  %14.7lg\n", itheta, jphi, chi[i][j] ) ;
    }
  }
  fclose(fp) ;
  /* NO phi-average
  for ( j = 0 ; j < nphi ; j++ ) chi[0][j] = 0. ;
  for ( i = 1 ; i < ntheta ; i++ ) {
    chimin = chi[i][0] ;
    chimax = chimin ;
    chiav = chimin ;
    for ( j = 1 ; j < nphi - 1 ; j++ ) {
      chiij = chi[i][j] ;
      chiav += chiij ; 
      if ( chiij < chimin ) chimin = chiij ;
      else if ( chiij > chimax ) chimax = chiij ; 
    }
    chiav /= (double) ( nphi - 1 ) ;
    chivar = chimax - chimin ;
    if ( chivar > EPSILON ) 
      for ( j = 0 ; j < nphi ; j++ ) 
        chi[i][j] = ( chi[i][j] - chiav ) / chivar ;
    else
      for ( j = 0 ; j < nphi ; j++ ) 
        chi[i][j] -=  chiav ;	    
  }
  end NO phi-average */
  chimin = chi[0][0] ;
  chimax = chimin ;
  for ( i = 0 ; i < ntheta ; i++ ) 
    for ( j = 0 ; j < nphi ; j++ ) {
      chiij = chi[i][j] ;
      if ( chiij < chimin ) chimin = chiij ;
      else if ( chiij > chimax ) chimax = chiij ; 
    }

  chivar = ( chimax - chimin ) / 45.999999999999 ;
  printf("chi min=%lf max=%lf var=%lf\n", chimin, chimax, chivar ) ;

  if ( chiminin != 1. || chimaxin != -1. ) { 
    if ( chiminin !=  1. ) chimin = chiminin ;
    if ( chimaxin != -1. ) chimax = chimaxin ;
    chivar = ( chimax - chimin ) / 45.999999999999 ;
    printf("REDEFINED: chi min=%lf max=%lf var=%lf\n", 
             chimin, chimax, chivar ) ;
  }
  for ( i = 0 ; i < ntheta ; i++ ) {
    for ( j = 0 ; j < nphi ; j++ ) {
      chi[i][j] = ( chi[i][j] - chimin ) / chivar ; /* normal */
      /* chi[i][j] = ( chimax - chi[i][j] ) / chivar ;  reverse */
      /* printf("%2d ", (int) chi[i][j] )	;  */
    }
    /* printf("\n") ; */
  }
/*
  npoints = POINTSPERDEG * 2 * ( thetamax + MARGIN ) + 1 ;
*/
  npoints = POINTSPERDEG * 2 * THETAFRAME + 1 ;

  printf("thetamax = %lf, npoints = %d\n", thetamax, npoints ) ;

  chilinelength = npoints ;
  chiline = ( int * ) calloc( chilinelength, sizeof( int ) ) ; 

  sprintf(pngfile,"%s_%.0lf_%.0lf.png",argv[1],1.e4*chimin,1.e4*chimax) ;
  fp = fopen("temp.xpm","w") ;
  fprintf(fp,"/* XPM */\n") ;
  fprintf(fp,"/* Black: %lf White: %lf */\n", chimin, chimax ) ;
  fprintf(fp,"static char * pict[]= {\n\"") ;
  fprintf(fp,"%d %d %d %d\",\n", npoints, npoints, 47, 1 ) ;
  for ( i = 0, j = 0 ; i < 15 ; i++, j++ ) 
    fprintf(fp,"\"%c c #%2s0000 \",", ( j<26 ? 'a'+j : 'A'+j-26 ), col16[i] ) ;
  fprintf(fp,"\n") ;
  for ( i = 0 ; i < 15 ; i++, j++ ) 
    fprintf(fp,"\"%c c #ff%2s00 \",", ( j<26 ? 'a'+j : 'A'+j-26 ), col16[i] ) ;
  fprintf(fp,"\n") ;
  for ( i = 0 ; i < 16 ; i++, j++ ) 
    fprintf(fp,"\"%c c #ffff%2s \",", ( j<26 ? 'a'+j : 'A'+j-26 ), col16[i] ) ;
  fprintf(fp,"\n\"U c None\"\n}\n") ;
  iy80  = ( THETAFRAME - 80 ) * POINTSPERDEG ;
  imy80 = ( THETAFRAME + 80 ) * POINTSPERDEG ;
  jx80  = ( THETAFRAME + 80 ) * POINTSPERDEG ;
  jmx80 = ( THETAFRAME - 80 ) * POINTSPERDEG ;
  for ( ipt = 0 ; ipt < npoints ; ipt++ ) {
    y = THETAFRAME - ( ipt / (double) POINTSPERDEG ) ;
    for ( jpt = 0 ; jpt < npoints ; jpt++ ) {
      x = ( jpt / (double) POINTSPERDEG ) - THETAFRAME ;
      r = sqrt( x*x + y*y )/90. ;
/*
      th = ( r < 1. ? asin( 2*r/(1.+r*r) ) : PI-asin( 2*r/(1.+r*r) ) ) ;
*/
      th = 2*atan(r) ;
      th *= 180./PI ; 
      if ( th < thetamax ) {
        ph = atan2( y, x ) ; 
        ph *= 180./PI ; 
        if ( ph < 0. ) ph += 360 ;
        chichi = interpol( th, ph, ntheta, dtheta, nphi, dphi, chi ) ;
	if ( chichi < 0. ) ichichi = 0 ;
	else if ( chichi > 45.999999999 ) ichichi = 45 ;
	else ichichi = (int) chichi ;
      }
      /* else if ( fabs( th - 90. ) < .65/POINTSPERDEG ) ichichi = 0 ; */
      else ichichi = 46 ;
      chiline[ jpt ] = ichichi ;
/*      printf("th = %lf, ph = %lf, chi = %lf\n", th, ph, chichi ) ; */
    }
    writechiline(fp, chilinelength, chiline ) ;
  }
  fprintf(fp,"} ;\n") ;
  fclose(fp) ;
  sprintf(command,"convert temp.xpm %s", pngfile ) ; system(command) ;

  return 0 ;
}

double interpol( double theta, double phi, int ntheta, double dtheta,
                 int nphi, double dphi, double **chi ) 
{
  int i1, j1, i2, j2 ;
  double th1, th2, ph1, ph2 ;
   /* find points (th1,ph1), (th1,ph2), (th2,ph1), (th2,ph2) */
/*
  printf("theta %lf, phi %lf, dtheta %lf, dphi %lf\n",
         theta, phi, dtheta, dphi ) ;
*/
  i1 = (int) ( theta / dtheta ) ;
  i2 = i1 + 1 ;
  j1 = (int) ( phi / dphi ) ;
  j2 = j1 + 1 ;  
  th2 = theta / dtheta - i1 ;
  th1 = 1. - th2 ;
  ph2 = phi / dphi - j1 ; 
  ph1 = 1. - ph2 ; 
/*
  printf("i1 %d, i2 %d, j1 %d, j2 %d, th1 %lf, th2 %lf, ph1 %lf, ph2 %lf\n",
          i1, i2, j1, j2, th1, th2, ph1, ph2 ) ;
*/
  return ( chi[i1][j1] * th1 * ph1 + chi[i1][j2] * th1 * ph2 
         + chi[i2][j1] * th2 * ph1 + chi[i2][j2] * th2 * ph2 ) ;
}

void writechiline(FILE *fp, int chilinelength, int *chiline ) 
{
  int i, j ;
  fprintf(fp,"\"") ;
  for ( i = 0 ; i < chilinelength ; i++ ) {
    j = chiline[i] ;
    fprintf(fp,"%c", ( j < 26 ? 'a'+j : 'A'+j-26 ) ) ;
  }
  fprintf(fp,"\",\n") ;
}

/* 16 colors 
void writechiline(FILE *fp, int chilinelength, int *chiline ) 
{
  static char cchi[17] =
    { '0','1','2','3','4','5','6','7','8','9','a','b','c','d','e','f','g'} ;

  int i ;
  fprintf(fp,"\"") ;
  for ( i =0 ; i < chilinelength ; i++ )
    fprintf(fp,"%c", cchi[ chiline[i] ]) ;
  fprintf(fp,"\",\n") ;
}
*/

/*
  fprintf(fp,"\"a c #000000 \",\n") ;
  fprintf(fp,"\"b c #110000 \",\n") ;
  fprintf(fp,"\"c c #220000 \",\n") ;
  fprintf(fp,"\"d c #330000 \",\n") ;
  fprintf(fp,"\"e c #440000 \",\n") ;
  fprintf(fp,"\"f c #550000 \",\n") ;
  fprintf(fp,"\"g c #660000 \",\n") ;
  fprintf(fp,"\"h c #770000 \",\n") ;
  fprintf(fp,"\"i c #880000 \",\n") ;
  fprintf(fp,"\"j c #990000 \",\n") ;
  fprintf(fp,"\"k c #aa0000 \",\n") ;
  fprintf(fp,"\"l c #bb0000 \",\n") ;
  fprintf(fp,"\"m c #cc0000 \",\n") ;
  fprintf(fp,"\"n c #dd0000 \",\n") ;
  fprintf(fp,"\"o c #ee0000 \",\n") ;
  fprintf(fp,"\"p c #ff0000 \",\n") ;
  fprintf(fp,"\"q c #ff1100 \",\n") ;
  fprintf(fp,"\"r c #ff2200 \",\n") ;
  fprintf(fp,"\"s c #ff3300 \",\n") ;
  fprintf(fp,"\"t c #ff4400 \",\n") ;
  fprintf(fp,"\"u c #ff5500 \",\n") ;
  fprintf(fp,"\"v c #ff6600 \",\n") ;
  fprintf(fp,"\"w c #ff7700 \",\n") ;
  fprintf(fp,"\"x c #ff8800 \",\n") ;
  fprintf(fp,"\"y c #ff9900 \",\n") ;
  fprintf(fp,"\"z c #ffaa00 \",\n") ;
  fprintf(fp,"\"A c #ffbb00 \",\n") ;
  fprintf(fp,"\"B c #ffcc00 \",\n") ;
  fprintf(fp,"\"C c #ffdd00 \",\n") ;
  fprintf(fp,"\"D c #ffee00 \",\n") ;
  fprintf(fp,"\"E c #ffff00 \",\n") ;
  fprintf(fp,"\"F c #ffff11 \",\n") ;
  fprintf(fp,"\"G c #ffff22 \",\n") ;
  fprintf(fp,"\"H c #ffff33 \",\n") ;
  fprintf(fp,"\"I c #ffff44 \",\n") ;
  fprintf(fp,"\"J c #ffff55 \",\n") ;
  fprintf(fp,"\"K c #ffff66 \",\n") ;
  fprintf(fp,"\"L c #ffff77 \",\n") ;
  fprintf(fp,"\"M c #ffff88 \",\n") ;
  fprintf(fp,"\"N c #ffff99 \",\n") ;
  fprintf(fp,"\"O c #ffffaa \",\n") ;
  fprintf(fp,"\"P c #ffffbb \",\n") ;
  fprintf(fp,"\"Q c #ffffcc \",\n") ;
  fprintf(fp,"\"R c #ffffdd \",\n") ;
  fprintf(fp,"\"S c #ffffee \",\n") ;
  fprintf(fp,"\"T c #ffffff \",\n") ;
*/
/*
  fprintf(fp,"%d %d %d %d\",\n", npoints, npoints, 17, 1 ) ;
  fprintf(fp,"\"0 c #000000\",\n") ;
  fprintf(fp,"\"1 c #330000\",\n") ;
  fprintf(fp,"\"2 c #660000\",\n") ;
  fprintf(fp,"\"3 c #990000\",\n") ;
  fprintf(fp,"\"4 c #cc0000\",\n") ;
  fprintf(fp,"\"5 c #ff0000\",\n") ;
  fprintf(fp,"\"6 c #ff3300\",\n") ;
  fprintf(fp,"\"7 c #ff6600\",\n") ;
  fprintf(fp,"\"8 c #ff9900\",\n") ;
  fprintf(fp,"\"9 c #ffcc00\",\n") ;
  fprintf(fp,"\"a c #ffff00\",\n") ;
  fprintf(fp,"\"b c #ffff33\",\n") ;
  fprintf(fp,"\"c c #ffff66\",\n") ;
  fprintf(fp,"\"d c #ffff99\",\n") ;
  fprintf(fp,"\"e c #ffffcc\",\n") ;
  fprintf(fp,"\"f c #ffffff\",\n") ;
  fprintf(fp,"\"g c None\"\n}\n") ;
  */
  
