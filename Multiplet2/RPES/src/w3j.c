#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "w3j.h"

#define LMAX  4
#define SDMAX 8
#define JDMAX (imax( 2*LMAX, SDMAX ))
#define iabs(i) ( (i) > 0 ? (i) : -(i) )
#define iodd(i) ( (i) - 2 * ( (i)/2 ) )
#define imax(i,j) ( (i) > (j) ? (i) : (j) )
#define imin(i,j) ( (i) < (j) ? (i) : (j) )


typedef struct {
    Jmind jm[5];
    double val;
} W3jord;

int w3jtabdim;
W3jord* w3jtab;


void iswap(int* i1, int* i2){ int idum = *i1; *i1 = *i2; *i2 = idum; }

int Jmind5less(Jmind* jmp1, Jmind* jmp2){
    int i;
    for (i=0; i<5; i++){
	if (jmp1[i] < jmp2[i]) return 1;
        else if (jmp2[i] < jmp1[i]) return 0;
    }
    return 0;
}

int Jmind5equal(Jmind* jmp1, Jmind* jmp2){
    int i;
    for (i=0; i<5; i++){
	if (jmp1[i] != jmp2[i]) return 0;
    }
    return 1;
}

int w3jtabindex(int jd1, int jd2, int jd3, int md1, int md2){
/* looks up (ordered) w3j in w3jtab, returns table index */
    Jmind jmp[5];
    int ilow, ihigh, imid;

    jmp[0] = (Jmind) jd1 ;
    jmp[1] = (Jmind) jd2 ; 
    jmp[2] = (Jmind) jd3 ; 
    jmp[3] = (Jmind) md1 ;
    jmp[4] = (Jmind) md2 ;
    ilow = imid = 0;
    ihigh = w3jtabdim - 1;

    while (ilow < ihigh){
	if (Jmind5less(jmp, w3jtab[imid].jm)) ihigh = imid - 1;
	else if (Jmind5less(w3jtab[imid].jm, jmp)) ilow = imid + 1;
	else break;
	imid = (ilow + ihigh)/2;
    }
    if (Jmind5equal(jmp, w3jtab[imid].jm))
	return imid;
    else
	return -1;   /* not found */
}

double w3jtabval(int jd1, int jd2, int jd3, int md1, int md2){
/* looks up (ordered) w3j in w3jtab, returns value */
    Jmind jmp[5];
    int ilow, ihigh, imid;

    jmp[0] = (Jmind) jd1 ;
    jmp[1] = (Jmind) jd2 ; 
    jmp[2] = (Jmind) jd3 ; 
    jmp[3] = (Jmind) md1 ;
    jmp[4] = (Jmind) md2 ;
    ilow = imid = 0;
    ihigh = w3jtabdim - 1;

    while (ilow < ihigh){
	if (Jmind5less(jmp, w3jtab[imid].jm)) ihigh = imid - 1;
	else if (Jmind5less(w3jtab[imid].jm, jmp)) ilow = imid + 1;
	else break;
	imid = (ilow + ihigh)/2;
    }
    if (Jmind5equal(jmp, w3jtab[imid].jm))
	return w3jtab[imid].val ;
    else {
	printf("error in w3jtabval: (jd1 ... md2) not found\n");
	exit(1);
    }
}

int w3jrestrict(int jd1, int jd2, int jd3, int md1, int md2, int md3){
/*  restrictions to odometer ordering of w3j symbols while writing w3jtab
    returns 0, 1  for this combination allowed, or not (=restricted), resp. */
    return md1 < 0 || md1 == 0 && md2 < 0 || jd1 == jd2 && iabs(md1) < iabs(md2)
	   || jd2 == jd3 && md2 < md3 || jd1 == jd3 && iabs(md1) < iabs(md3) ;
}

int w3jdefined(int jd1, int jd2, int jd3, int md1, int md2, int md3){
/* returns 0, 1 for 3j-symbol undefined, defined, respectively */
    if ( jd1 < 0 || iabs(md1) > jd1 || iodd( jd1 + md1 ) || 
	 jd2 < 0 || iabs(md2) > jd2 || iodd( jd2 + md2 ) ||
         jd3 < 0 || iabs(md3) > jd3 || iodd( jd3 + md3 )   )
       return 0;
    else 
       return 1;
}

int w3jnonzero(int jd1, int jd2, int jd3, int md1, int md2, int md3){
/* returns 0, 1 for 3j-symbol zero or nonzero, respectively */
    if ( jd1 + jd2 < jd3 || jd2 + jd3 < jd1  || jd3 + jd1 < jd2 ||
	      iodd( jd1 + jd2 + jd3 ) ||  iodd( jd1 - jd2 - md3 ) || 
	      md1 + md2 + md3 != 0 || iodd( (jd1+jd2+jd3)/2 ) && 
              ( md1 == 0 && md2 == 0 || jd1 == jd2 && md1 == md2 ||
	        jd1 == jd3 && md1 == md3 || jd2 == jd3 && md2 == md3 ) )
       return 0;
    else
       return 1;
}

int w3jorder(int* jd1, int* jd2, int* jd3, int* md1, int* md2, int* md3){
/* brings (jd1,jd2,jd3,md1,md2,md3) in odometer order and returns sign change */
    int mc = 0;
    if ( *jd1 < *jd2) { iswap(jd1,jd2); iswap(md1,md2); mc++; }
    if ( *jd2 < *jd3) { iswap(jd2,jd3); iswap(md2,md3); mc++; }
    if ( *jd1 < *jd2) { iswap(jd1,jd2); iswap(md1,md2); mc++; }
    if ( *jd1 == *jd2 && iabs(*md1) < iabs(*md2) ) { iswap(md1,md2); mc++; }
    if ( *jd1 == *jd3 && iabs(*md1) < iabs(*md3) ) { iswap(md1,md3); mc++; }
    if ( *md1 < 0 || *md1 == 0 && *md2 < 0 ){
	*md1 *= -1; *md2 *= -1; *md3 *= -1; mc++; }
    if ( *jd2 == *jd3 && *md2 < *md3 ) { iswap(md2,md3); mc++; }
    if ( iodd(mc) && iodd( ( *jd1 + *jd2 + *jd3 )/2 ) )
       return -1;
    else
       return 1;
}					

double w3jcalc(int jd1, int jd2, int jd3, int md1, int md2, int md3){
/* calculates w3j-symbol using closed factorial formula and log(n!)'s */
    int iabcp, iabc, icab, ibca, iapd, iamd, ibpe, ibme, icpf, icmf ;
    int nzmic2, nzmic3, nzmi, nzmx, nz, nzt1, nzt2, nzt3, nzt4, nzt5 ;
    int i, maxnfac ;
    double sqfclg, s1, termlg, ssterm, val;
    double *faclog;

    if (!w3jdefined(jd1, jd2, jd3, md1, md2, md3)){
       printf(" !w3jdefined in w3jcalc\n"); 
       exit(1);
    }
    else if (!w3jnonzero(jd1, jd2, jd3, md1, md2, md3))
       return 0.e0;
    else if ( jd1 == 0 || jd2 == 0 ) 
	return (iodd((jd3-md3)/2) ? -1. : 1.) / sqrt( (double) (jd3+1) ) ;
    else if ( jd3 == 0 )
	return (iodd((jd1-md1)/2) ? -1. : 1.) / sqrt( (double) (jd2+1) ) ;
    else {
	iabcp = ( jd1 + jd2 + jd3 )/2 ;
	iabc  = iabcp - jd3 ;
	icab  = iabcp - jd2 ;
	ibca  = iabcp - jd1 ;
	iapd  = ( jd1 + md1 )/2 ;
	iamd  = iapd - md1 ;
	ibpe  = ( jd2 + md2 )/2 ;
	ibme  = ibpe - md2 ;
	icpf  = ( jd3 - md3 )/2 ;
	icmf  = icpf + md3 ;
	nzmic2 = (jd2 - jd3 - md1 )/2 ;
	nzmic3 = (jd1 - jd3 + md2 )/2 ;
	nzmi = imax(imax(0,nzmic2),nzmic3) ;
	nzmx = imin(imin(iabc,iamd),ibpe) ; 
	s1 = iodd(nzmi) ? -1. : 1. ;   
	if ( nzmx < nzmi )
	    return 0.; 
	maxnfac =  iabcp + 2 ;
	faclog = ( double * ) calloc( maxnfac, sizeof(double) ) ;  
	faclog[0] = faclog[1] = 0. ;
	for ( i = 2 ; i < maxnfac  ; i++ )
	    faclog[i] = faclog[i-1] + log( (double) i ) ;
	sqfclg= 0.5 * (- faclog[iabcp+1] 
		       + faclog[iabc] + faclog[icab] + faclog[ibca]
		       + faclog[iapd] + faclog[iamd] + faclog[ibpe]
		       + faclog[ibme] + faclog[icpf] + faclog[icmf] ) ;
	val = 0. ;
	for ( nz = nzmi ; nz <= nzmx ; nz++ ) {
	    nzt1 = iabc - nz ;
	    nzt2 = iamd - nz ;
	    nzt3 = ibpe - nz ;
	    nzt4 = nz - nzmic2 ;
	    nzt5 = nz - nzmic3 ;
	    termlg = sqfclg - faclog[nz] - faclog[nzt1] - faclog[nzt2]
		          - faclog[nzt3] - faclog[nzt4] - faclog[nzt5] ;
	    ssterm = s1 * exp(termlg) ;
	    val = val + ssterm ;
	    s1 = -s1 ;
	}
	val *= (iodd((jd1-jd2-md3)/2) ? -1. : 1.) ;
	free(faclog) ;
	return val ;
    }
}

void w3jtabmake(){
/* 
   defines and allocates table (W3jord*) w3jtab, 
   and its dimension (int) w3jtabdim 
*/
    int jd1,jd2,jd3,md1,md2,md3,count ;
/*
    count w3jtabdim
*/
    count = 0 ;
    for (jd1=0;jd1<=JDMAX;jd1++){
     for (jd2=0;jd2<=jd1;jd2++){
      for (jd3=iabs(jd1-jd2);jd3<=jd2;jd3+=2){
       if ( imax(imax(jd1,jd2),jd3) <= SDMAX ||
	    !iodd(jd1) && !iodd(jd2) ) {
        for (md1=jd1-2*(jd1/2);md1<=jd1;md1+=2){
         for (md2=-jd2;md2<=jd2;md2+=2){
	     md3 = - md1 - md2;
	     if (iabs(md3) <= jd3) {
		 if ( !w3jrestrict(jd1,jd2,jd3,md1,md2,md3) 
		      && w3jnonzero(jd1,jd2,jd3,md1,md2,md3) )
		      ++count ;
	     }
         }   
        }
       }	
      }
     }	
    }
    w3jtabdim = count ;

/* allocate w3jtab */

    w3jtab = ( W3jord * ) calloc( w3jtabdim, sizeof(W3jord) ) ;

/* calculate w3j and write into w3jtab */

    count = 0 ;
    for (jd1=0;jd1<=JDMAX;jd1++){
     for (jd2=0;jd2<=jd1;jd2++){
      for (jd3=iabs(jd1-jd2);jd3<=jd2;jd3+=2){
       if ( imax(imax(jd1,jd2),jd3) <= SDMAX ||
	    !iodd(jd1) && !iodd(jd2) ) {
        for (md1=jd1-2*(jd1/2);md1<=jd1;md1+=2){
         for (md2=-jd2;md2<=jd2;md2+=2){
	     md3 = - md1 - md2;
	     if (iabs(md3) <= jd3) {
		 if ( !w3jrestrict(jd1,jd2,jd3,md1,md2,md3) 
		      && w3jnonzero(jd1,jd2,jd3,md1,md2,md3) ) {
		     w3jtab[count].jm[0] = (Jmind) jd1;
		     w3jtab[count].jm[1] = (Jmind) jd2;
		     w3jtab[count].jm[2] = (Jmind) jd3;
		     w3jtab[count].jm[3] = (Jmind) md1;
		     w3jtab[count].jm[4] = (Jmind) md2;
		     w3jtab[count].val = w3jcalc(jd1,jd2,jd3,md1,md2,md3) ;
/*
		     if ( jd1 < 9 && jd2 < 9 && jd3 < 9 )
			 printf("%5d: ( %2d %2d %2d |%3d %3d %3d ) = %.8lf\n", 
			    count,jd1,jd2,jd3,md1,md2,md3,w3jtab[count].val) ; 
*/
		     ++count ;
	       }
	     }
         }   
        }
       }	
      }
     }	
    }
    printf(" w3jtabmake: LMAX = %d, SDMAX = %d, JDMAX = %d,",LMAX,SDMAX,JDMAX);
    printf(" count = %d, w3jtabdim = %d\n", count, w3jtabdim);
}

double w3j(int jd1, int jd2, int jd3, int md1, int md2, int md3){
/*
  looks up value of w3j using w3jtabval if in range of table,
  otherwise calculates it using w3jcalc (see above)
*/
    int jd123max, sign ;

    if (!w3jdefined(jd1, jd2, jd3, md1, md2, md3)){
       printf(" !w3jdefined in w3j\n"); 
       exit(1);
    }
    else if (!w3jnonzero(jd1, jd2, jd3, md1, md2, md3))
       return 0.e0;
    else {
	jd123max = imax(imax(jd1,jd2),jd3) ;
	if ( jd123max <= SDMAX || 
	     jd123max <= JDMAX && !iodd(jd1) && !iodd(jd2) ){
	    sign = w3jorder(&jd1, &jd2, &jd3, &md1, &md2, &md3);
	    return sign * w3jtabval(jd1,jd2,jd3,md1,md2);
	}
	else
	    return w3jcalc(jd1, jd2, jd3, md1, md2, md3);
    }
}

double clebsch(int jd1, int md1, int jd2, int md2, int jd3, int md3){
/* 
   calculates clebsch-gordon coefficient from wigner 3j symbol, using w3j
   phase convention according to Cowan
*/
    return (iodd((jd1-jd2+md3)/2) ? -1. : 1.) 
	* sqrt( (double) (jd3+1) )
	* w3j(jd1, jd2, jd3, md1, md2, -md3) ;
}


