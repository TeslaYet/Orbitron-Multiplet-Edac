#include <stdio.h>
#include "globals.h"

void fockcalcs( struct Fock *state ) {
  int k, sign ;
  state -> s = ZERO ;
  sign = 1 ;
  for ( k = NBITS - 1 ; k >= 0 ; k-- ) {
    if ( fockocc( state, k ) == 1 ) 
      sign *= -1 ;
    if ( sign == -1 ) 
      state -> s |= ONE << k  ;
  }
}
void fockinitdual( struct Fock *state, Dualrep nin ) { 
  state -> n = nin ; fockcalcs( state ) ; 
}
void fockinit( struct Fock *state, int norb, int *occorb ) {
  int i ;
  state -> n = ZERO ;
  for ( i = 0 ; i < norb ; i++ )
    state -> n |= ONE << occorb[i] ;
  fockcalcs( state ) ;
}
int fockocc( struct Fock *state, int iorb ) {
  return (int) ( ( state -> n >> iorb ) & ONE ) ;
}
int fockcre( struct Fock *state, int iorb ) {
  Dualrep maskend, sendneg ;
  if ( ( state -> n >> iorb ) & ONE == ONE )  
    return 0 ;
  else {
    state -> n |= ONE << iorb ;
    maskend = ~ZERO << ( iorb + 1 ) ;
    sendneg = ~( state -> s | maskend ) ;
    state -> s = state -> s & maskend | sendneg ;
    return ( ( state -> s >> (iorb+1) ) & ONE == ONE ) ? -1 : 1  ; 
  }
}
int fockdes( struct Fock *state, int iorb ) {
  Dualrep maskend, sendneg ;
  if ( ( int ) ( ( state -> n >> iorb ) & ONE ) == 0 )  
    return 0 ;
  else {
    state -> n &= ~( ONE << iorb ) ;
    maskend = ~ZERO << ( iorb + 1 ) ;
    sendneg = ~( state -> s | maskend ) ;
    state -> s = state -> s & maskend | sendneg ;
    return ( ( state -> s >> (iorb+1) ) & ONE == ONE ) ? -1 : 1  ; 
  }
}
void fockdisplay( struct Fock *state ) {
  int iorb ;
  printf("|") ;
  for ( iorb = NBITS-1 ; iorb >= 0 ; iorb-- )
    printf("%d", fockocc( state, iorb ) ) ;
  printf(">") ;
}
void fockdisplays ( struct Fock *state ) {
  int iorb ;
  printf("|") ;  
  for ( iorb = NBITS-1 ; iorb >= 0 ; iorb-- )
    printf("%d", (int) ( ( state -> s >> iorb ) & ONE ) ) ;
  printf(">") ;
}  

void focklinearsort( int nstates, struct Fock *state ) 
{
  int i, imin, j ;
  Dualrep min ;
  struct Fock dummystate, *pstate ;

  for ( i = 0 ;  i < nstates - 1 ; i++ ) {
    pstate = state + i ;
    imin = i ;
    min = pstate -> n ;
    for ( j = i + 1 ; j < nstates ; j++ )
      if ( (++pstate) -> n < min ) { 
	imin = j ;
	min = pstate -> n ;
      }
    dummystate = state[i] ;
    state[i] = state[imin] ;
    state[imin] = dummystate ;
  }
}

/*****************************************************************************/
/* findstate finds the array position k of some state "*statek" in the       */
/*           array "state0" of dimension nstates. Uses binary search.        */
/*           &k  on exit such that  state0[k-1].n < statek->n <= state0[k].n */
/*           Returns  1  if  *statek  is found  0  otherwise                 */
/*****************************************************************************/
int findstate( int nstates, struct Fock *state0, struct Fock *statek, int *k ) 
{
  int khigh, klow, kmid ;
  Dualrep nk ;
  
  nk = statek -> n ;
  khigh = nstates ;
  klow = -1 ;
  while ( khigh > klow + 1 ) {
    kmid = ( klow + khigh ) / 2 ;
    if (  nk  <=  state0[kmid].n )
      khigh = kmid ;
    else                                    
      klow = kmid ;
  }
  *k = khigh ;
  if ( nk == state0[khigh].n )
    return 1 ;
  else
    return 0 ;
}
/* we always have klow < khigh && state0[klow].n < nj <= state0[khigh].n  */
/* thus at the end, when khigh = klow + 1, the right state is No khigh    */
