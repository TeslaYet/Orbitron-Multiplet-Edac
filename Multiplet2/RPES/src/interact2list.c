#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "globals.h"
#include "w3j.h"

int sporbindex( int *sorb1sh, int *lsh, int ish, int m, int sigma )
{
  return sorb1sh[ish] + 2 * ( lsh[ish] + m ) + ( sigma > 0 ? 1 : 0 ) ;
}

double cklmlm( int k, int l1, int m1, int l2, int m2 ) {
  if ( !w3jdefined( 2*l1, 2*k, 2*l2, -2*m1, 2*(m1-m2), 2*m2 ) ||
       !w3jnonzero( 2*l1, 2*k, 2*l2, -2*m1, 2*(m1-m2), 2*m2 ) )
    return 0. ;
  else 
    return ( m1 % 2 == 0 ? 1 : -1 )              /* = (-)^m == (-)^{-m} */
      * sqrt( ( 2 * l1 + 1 ) * ( 2 * l2 + 1 ) )  /* = [l1,l2]^{1/2} */
      * w3j( 2*l1, 2*k, 2*l2, 0, 0, 0 )
      * w3j( 2*l1, 2*k, 2*l2, -2*m1, 2*(m1-m2), 2*m2 ) ; 
}


int coulomb2list( int nshells, int *lsh, int *sorb1sh, 
		  struct CIlistitem *pcilistitem0, 
		  struct O2plistitem **ppo2plistitem0 )
{
  int ish1, ish2, ish3, ish4, l1, l2, l3, l4, m1, m2, m3, m4, s1, s2, s3, s4 ;
  int i1, i2, i3, i4, nk, k, counter = 0 ;
  double *r4312, v0 ;
  struct CIlistitem *ip = pcilistitem0 ;

  *ppo2plistitem0 = 0 ;

  while( ip != 0 ) {
    ish1 = ip->ish1 ; ish2 = ip->ish2 ; ish3 = ip->ish3 ; ish4 = ip->ish4 ; 
    l1 = lsh[ish1] ; l2 = lsh[ish2] ; l3 = lsh[ish3] ; l4 = lsh[ish4] ;
    nk = ip -> nk ; /* nk=min(l4+l1,l2+l3)+1. ck(4,1)ck(2,3)!=0 for k=0..nk-1*/
    r4312 = ip -> rmx ; /* r4312[k] defined for k=0..nk-1 */

    for ( m1 = -l1 ; m1 <= l1 ; m1++ ) 
      for ( m2 = -l2 ; m2 <= l2 ; m2++ )
	for ( m3 = -l3 ; m3 <= l3 ; m3++ ) { 
	  m4 = m1 + m2 - m3 ;
	  if ( m4 >= -l4 && m4 <= l4 ) 
	    for ( s1 = -1 ; s1 <= 1 ; s1 += 2 )
	      for ( s2 = -1 ; s2 <= 1 ; s2 += 2 ) {
		s4 = s1 ; 
		s3 = s2 ;
		i1 = sporbindex( sorb1sh, lsh, ish1, m1, s1 ) ;
		i2 = sporbindex( sorb1sh, lsh, ish2, m2, s2 ) ;
		i3 = sporbindex( sorb1sh, lsh, ish3, m3, s3 ) ;
		i4 = sporbindex( sorb1sh, lsh, ish4, m4, s4 ) ;
		if ( i1 < i2 && i3 != i4 ) { 
		  v0 = 0. ;
		  for ( k = 0 ; k < nk ; k++ )
		    v0 += cklmlm(k,l4,m4,l1,m1) * cklmlm(k,l2,m2,l3,m3) 
		          * r4312[k] ;
		  if ( v0 != 0. ) {
		    o2plistitemadd( ppo2plistitem0, i1, i2, i3, i4, v0 ) ;
		    counter++ ;
		  }
		}

	      }
	}
    ip = ip -> next ;
  }
  return counter ;
}

int vai2list( int nshells, int *lsh, int *sorb1sh, 
              struct O2plistitem **ppo2plistitem0 )
{
  int ish1, ish2, ish3, ish4, l1, l2, l3, l4, m1, m2, m3, m4, s1, s2, s3, s4 ;
  int n1, n2, i1, i2, i3, i4, nk, k, counter = 0 ;
  double *r4312, v0 ;
  *ppo2plistitem0 = 0 ;

/* < ish4, ish3 | V| ish1, ish2 > c4+ c3+ c2 c1. Here <0,2-end|V|1,1>  */

  ish1 = ish2 = 1 ; ish4 = 0 ;   
  l1 = lsh[ish1] ; l2 = lsh[ish2] ; l4 = lsh[ish4] ; 
  n1 = l4 + l1 ; 

  for ( ish3 = 2 ; ish3 < nshells ; ish3++ ) {
    l3 = lsh[ish3] ; n2 = l2 + l3 ;
    /* nk=min(l4+l1,l2+l3)+1. ck(4,1)ck(2,3)!=0 for k=0..nk-1*/
    nk = ( n1 < n2 ? n1 : n2 ) + 1 ;
    r4312 = ( double * ) calloc( nk, sizeof( double ) ) ; 
    /* read or calculate r4312[k] == rk(4,3;1,2)  k=0..nk-1 */
    printf("Enter Rk(%d,%d;%d,%d) for k=0..%d: ", l4, l3, l1, l2, nk-1 ) ;
    for ( k = 0 ; k < nk ; k++ ) scanf("%lf", &r4312[k] ) ;

    for ( m1 = -l1 ; m1 <= l1 ; m1++ ) 
      for ( m2 = -l2 ; m2 <= l2 ; m2++ )
	for ( m3 = -l3 ; m3 <= l3 ; m3++ ) { 
	  m4 = m1 + m2 - m3 ;
	  if ( m4 >= -l4 && m4 <= l4 ) 
	    for ( s1 = -1 ; s1 <= 1 ; s1 += 2 )
	      for ( s2 = -1 ; s2 <= 1 ; s2 += 2 ) {
		s4 = s1 ; 
		s3 = s2 ;
		i1 = sporbindex( sorb1sh, lsh, ish1, m1, s1 ) ;
		i2 = sporbindex( sorb1sh, lsh, ish2, m2, s2 ) ;
		i3 = sporbindex( sorb1sh, lsh, ish3, m3, s3 ) ;
		i4 = sporbindex( sorb1sh, lsh, ish4, m4, s4 ) ;
		if ( i3 != i4 ) { 
		  v0 = 0. ;
		  for ( k = 0 ; k < nk ; k++ )
		    v0 += cklmlm(k,l4,m4,l1,m1) * cklmlm(k,l2,m2,l3,m3) 
		          * r4312[k] ;
		  if ( v0 != 0. ) {
                    printf("<%d%2d%2d,%d%2d%2d|%d%2d%2d,%d%2d%2d> = %11.6lf\n",  
                              l4,m4,s4,l3,m3,s3,l1,m1,s1,l2,m2,s2, v0 ) ;
		    o2plistitemadd( ppo2plistitem0, i1, i2, i3, i4, v0 ) ;
		    counter++ ;
		  }
		}

	      }
	}
        free( r4312 ) ;
  }
  return counter ;
}

int spinorbit2list( int nshells, int *lsh, int *sorb1sh, double *ksish, 
		    struct O1plistitem **ppo1plistitem0 ) 
{
  int ish, l0, m, sigma, i1, i2, counter = 0 ;
  double ksi0, v0 ;

  for ( ish = 0 ; ish < nshells ; ish++ ) {
    l0   = lsh[ish] ; 
    ksi0 = ksish[ish] ;
    if ( l0 != 0 && ksi0 != 0. ) {
      /* lz sz term (diagonal) */
      for ( m = -l0 ; m <= l0 ; m++ )
	if ( m != 0 ) {
	  v0 = 0.5 * m * ksi0 ;
	  for ( sigma = -1 ; sigma <= 1 ; sigma += 2 ) {
	    i1 = sporbindex( sorb1sh, lsh, ish, m, sigma ) ; 
	    o1plistitemadd( ppo1plistitem0, i1, i1, sigma * v0 ) ;
            /*	printf("%2d %2d %20.12e\n", i1, i1, sigma * v0 ) ; */
	    counter++ ;
	  }
	}
      /* 0.5 ( l+ s- + l- s+ ) term (off-diagonal) */
      for ( m = -l0 ; m < l0 ; m++ ) {
	v0 = 0.5 * sqrt( ( l0 - m ) * ( l0 + m + 1 ) ) * ksi0 ;
	i1 = sporbindex( sorb1sh, lsh, ish, m,    1 ) ; 
	i2 = sporbindex( sorb1sh, lsh, ish, m+1, -1 ) ; 
	o1plistitemadd( ppo1plistitem0, i1, i2, v0 ) ;
	o1plistitemadd( ppo1plistitem0, i2, i1, v0 ) ;
	/*  printf("%2d %2d %20.12e (l+s-) and ", i1, i2, v0 ) ;
	    printf("%2d %2d %20.12e (l-s+)\n",    i2, i1, v0 ) ;  */
	counter += 2 ;
      }
    }
  }
  return counter ;
}

/* calculates  magnetic field operator = -hz * jz */
int magneticfield2list( int nshells, int *lsh, int *sorb1sh, double hz,
			struct O1plistitem **ppo1plistitem0 ) 
{
  int ish, l0, m, sigma, i1, counter = 0 ;

  for ( ish = 0 ; ish < nshells ; ish++ ) {
    l0   = lsh[ish] ; 
    for ( m = -l0 ; m <= l0 ; m++ )
      for ( sigma = -1 ; sigma <= 1 ; sigma += 2 ) {
	i1 = sporbindex( sorb1sh, lsh, ish, m, sigma ) ; 
/* case of magnetic field operator = -hz * jz */
/*
	o1plistitemadd( ppo1plistitem0, i1, i1, -hz * ( m + 0.5 * sigma ) ) ;
        printf("%2d %2d %20.12e\n", i1, i1, -hz * ( m + 0.5 * sigma ) );
*/ 
/* case of magnetic field operator = -hz * 2 * sz */
	o1plistitemadd( ppo1plistitem0, i1, i1, -hz * sigma ) ;
	counter++ ;
      }
  }
  return counter ;
}

/* applies exchange field  -2 * hz * sz   to shell ish */
int exchangefield2list( int nshells, int *lsh, int *sorb1sh, int ish,
			double hz, struct O1plistitem **ppo1plistitem0 ) 
{
  int l0, m, sigma, i1, counter = 0 ;
  
  if ( ish < nshells ) {
    l0  = lsh[ish] ;
    for ( m = -l0 ; m <= l0 ; m++ )
      for ( sigma = -1 ; sigma <= 1 ; sigma += 2 ) {
	i1 = sporbindex( sorb1sh, lsh, ish, m, sigma ) ; 
	/* case of magnetic field operator = -hz * 2 * sz */
	o1plistitemadd( ppo1plistitem0, i1, i1, -hz * sigma ) ;
	counter++ ;
      }
  }
  else {printf("exchangefield2list: ish >= nshells\n");exit(1);}
  return counter ;
}
/* calculates magnetic field operator = -h.( l+2s) with h=(hx,0,hz) */
int magnfieldxz2list( int nshells, int *lsh, int *sorb1sh, double h, 
		      double theta, struct O1plistitem **ppo1plistitem0 ) 
{
  int ish, l0, m, sigma, i1, i2, counter = 0 ;
  double hx, hz, v0 ;

  hz = h * cos( theta ) ;  
  hx = h * sin( theta ) ;
  for ( ish = 0 ; ish < nshells ; ish++ ) {
    l0   = lsh[ish] ; 
    for ( m = -l0 ; m <= l0 ; m++ ) {
      for ( sigma = -1 ; sigma <= 1 ; sigma += 2 ) {
	i1 = sporbindex( sorb1sh, lsh, ish, m, sigma ) ; 
/* -hz * ( lz + 2*sz ) */
	o1plistitemadd( ppo1plistitem0, i1, i1, -hz * ( m + sigma ) ) ;
	counter++ ;
	/*        printf("%2d %2d %20.12e\n", i1, i1, -hz * ( m + sigma ) ); */
      }
    }
/* -hx * ( s+ + s- ) */
    for ( m = -l0 ; m <= l0 ; m++ ) {
      i1 = sporbindex( sorb1sh, lsh, ish, m,  1 ) ; 
      i2 = sporbindex( sorb1sh, lsh, ish, m, -1 ) ; 
      o1plistitemadd( ppo1plistitem0, i1, i2, -hx ) ;
      o1plistitemadd( ppo1plistitem0, i2, i1, -hx ) ;
      /*      printf("%2d %2d %20.12e\n", i2, i1, -hx ); */
      counter += 2 ;
    }
/* -hx * ( l+ + l- ) / 2 */
    for ( m = -l0 ; m < l0 ; m++ ) {
      v0 = -0.5 * sqrt( ( l0 - m ) * ( l0 + m + 1 ) ) * hx ;
      for ( sigma = -1 ; sigma <= 1 ; sigma += 2 ) {
	i1 = sporbindex( sorb1sh, lsh, ish, m,   sigma ) ; 
	i2 = sporbindex( sorb1sh, lsh, ish, m+1, sigma ) ; 
	o1plistitemadd( ppo1plistitem0, i1, i2, v0 ) ;
	o1plistitemadd( ppo1plistitem0, i2, i1, v0 ) ;
      }
      counter += 2 ;
    }	
  }
  return counter ;
}

/* exchange field operator = -h.(2s) with h=(hx,0,hz)  on shell No ish */
int exchfieldxz2list( int nshells, int *lsh, int *sorb1sh, int ish, double h, 
		      double theta, struct O1plistitem **ppo1plistitem0 ) 
{
  int l0, m, sigma, i1, i2, counter = 0 ;
  double hx, hz, v0 ;

  if ( ish >= nshells ) {
    printf("exchfieldxz2list: ish=%d >= nshells=%d. Stop.\n",ish,nshells);
    exit(1) ;
  }
  hz = h * cos( theta ) ;  
  hx = h * sin( theta ) ;

    l0   = lsh[ish] ; 
    for ( m = -l0 ; m <= l0 ; m++ ) {
      for ( sigma = -1 ; sigma <= 1 ; sigma += 2 ) {
	i1 = sporbindex( sorb1sh, lsh, ish, m, sigma ) ; 
/* -hz * ( 2*sz ) */
	o1plistitemadd( ppo1plistitem0, i1, i1, -hz * sigma ) ;
	counter++ ;
	/*        printf("%2d %2d %20.12e\n", i1, i1, -hz * sigma ); */
      }
    }
/* -hx * ( s+ + s- ) */
    for ( m = -l0 ; m <= l0 ; m++ ) {
      i1 = sporbindex( sorb1sh, lsh, ish, m,  1 ) ; 
      i2 = sporbindex( sorb1sh, lsh, ish, m, -1 ) ; 
      o1plistitemadd( ppo1plistitem0, i1, i2, -hx ) ;
      o1plistitemadd( ppo1plistitem0, i2, i1, -hx ) ;
      /*      printf("%2d %2d %20.12e\n", i2, i1, -hx ); */
      counter += 2 ;
    }

  return counter ;
}

/* defines level of 1st 2 shells (shells No 0,1) */
int esh012list( int nshells, int *lsh, int *sorb1sh, double esh0, double esh1,
               struct O1plistitem **ppo1plistitem0 ) 
{
  int i1, counter = 0 ;
  for ( i1 = 0 ; i1 < sorb1sh[1] ; i1++ ) {
    o1plistitemadd( ppo1plistitem0, i1, i1, esh0 ) ;
    /*    printf("%2d %2d %20.12e\n", i1, i1, esh0 ) ;*/
    counter++ ;
  }
  for ( i1 = sorb1sh[1] ; i1 < sorb1sh[2] ; i1++ ) {
    o1plistitemadd( ppo1plistitem0, i1, i1, esh1 ) ;
    /*    printf("%2d %2d %20.12e\n", i1, i1, esh1 ) ;*/
    counter++ ;
  }
  return counter ;
}

/* defines level of 1st shell (shell No 0) */
int esh02list( int nshells, int *lsh, int *sorb1sh, double esh0,
               struct O1plistitem **ppo1plistitem0 ) 
{
  int i1, counter = 0 ;
  for ( i1 = 0 ; i1 < sorb1sh[1] ; i1++ ) {
    o1plistitemadd( ppo1plistitem0, i1, i1, esh0 ) ;
    printf("%2d %2d %20.12e\n", i1, i1, esh0 ) ;
    counter++ ;
  }
  return counter ;
}
