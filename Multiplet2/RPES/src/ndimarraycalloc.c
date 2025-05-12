#include <stdlib.h>
#include <complex.h>
double *calloc1double(int n1){
  double *p;
  p = ( double * ) calloc( n1, sizeof( double ) );
  return p ;
}

void free1double(double *p){
  free( p ) ;
}

double **calloc2double(int n1, int n2){
  int i1;
  double **p;
  p = ( double ** ) calloc( n1, sizeof( double * ) );
  for ( i1 = 0 ; i1 < n1 ; i1++ ) {
    p[i1] = ( double * ) calloc( n2, sizeof( double ) ); 
  }
  return p ;
}

void free2double(int n1, double **p){
  int i1 ;  
  for ( i1 = 0 ; i1 < n1 ; i1++ ) {
    free( p[i1] ) ;
  }
  free( p ) ;
}

double ***calloc3double(int n1, int n2, int n3){
  int i1, i2;
  double ***p;
  p = ( double *** ) calloc( n1, sizeof( double ** ) );
  for ( i1 = 0 ; i1 < n1 ; i1++ ) {
    p[i1] = ( double ** ) calloc( n2, sizeof( double * ) ); 
    for ( i2 = 0 ; i2 < n2 ; i2++ ) {
      p[i1][i2] = ( double * ) calloc( n3, sizeof( double ) );
    }
  }
  return p ;
}

void free3double(int n1, int n2, double ***p){
  int i1, i2;  
  for ( i1 = 0 ; i1 < n1 ; i1++ ) {
    for ( i2 = 0 ; i2 < n2 ; i2++ ) {
      free( p[i1][i2] ) ;
    }
    free( p[i1] ) ;
  }
  free( p ) ;
}

double ****calloc4double(int n1, int n2, int n3, int n4){
  int i1, i2, i3;
  double ****p;
  p = ( double **** ) calloc( n1, sizeof( double *** ) );
  for ( i1 = 0 ; i1 < n1 ; i1++ ) {
    p[i1] = ( double *** ) calloc( n2, sizeof( double ** ) ); 
    for ( i2 = 0 ; i2 < n2 ; i2++ ) {
      p[i1][i2] = ( double ** ) calloc( n3, sizeof( double * ) );
      for ( i3 = 0 ; i3 < n3 ; i3++ ) {
	p[i1][i2][i3] = ( double * ) calloc( n4, sizeof( double ) );
      }	
    }
  }
  return p ;
}

void free4double(int n1, int n2, int n3, double ****p){
  int i1, i2, i3;  
  for ( i1 = 0 ; i1 < n1 ; i1++ ) {
    for ( i2 = 0 ; i2 < n2 ; i2++ ) {
      for ( i3 = 0 ; i3 < n3 ; i3++ ) {
	free( p[i1][i2][i3] ) ;
      }
      free( p[i1][i2] ) ;
    }
    free( p[i1] ) ;
  }
  free( p ) ;
}

complex ****calloc4cplx(int n1, int n2, int n3, int n4){
  int i1, i2, i3;
  complex ****p;
  p = ( complex **** ) calloc( n1, sizeof( complex *** ) );
  for ( i1 = 0 ; i1 < n1 ; i1++ ) {
    p[i1] = ( complex *** ) calloc( n2, sizeof( complex ** ) ); 
    for ( i2 = 0 ; i2 < n2 ; i2++ ) {
      p[i1][i2] = ( complex ** ) calloc( n3, sizeof( complex * ) );
      for ( i3 = 0 ; i3 < n3 ; i3++ ) {
	p[i1][i2][i3] = ( complex * ) calloc( n4, sizeof( complex ) );
      }	
    }
  }
  return p ;
}

void free4cplx(int n1, int n2, int n3, complex ****p){
  int i1, i2, i3;  
  for ( i1 = 0 ; i1 < n1 ; i1++ ) {
    for ( i2 = 0 ; i2 < n2 ; i2++ ) {
      for ( i3 = 0 ; i3 < n3 ; i3++ ) {
	free( p[i1][i2][i3] ) ;
      }
      free( p[i1][i2] ) ;
    }
    free( p[i1] ) ;
  }
  free( p ) ;
}

float *calloc1float(int n1){
  float *p;
  p = ( float * ) calloc( n1, sizeof( float ) );
  return p ;
}

void free1float(float *p){
  free( p ) ;
}

float **calloc2float(int n1, int n2){
  int i1;
  float **p;
  p = ( float ** ) calloc( n1, sizeof( float * ) );
  for ( i1 = 0 ; i1 < n1 ; i1++ ) {
    p[i1] = ( float * ) calloc( n2, sizeof( float ) ); 
  }
  return p ;
}

void free2float(int n1, float **p){
  int i1 ;  
  for ( i1 = 0 ; i1 < n1 ; i1++ ) {
    free( p[i1] ) ;
  }
  free( p ) ;
}

float ***calloc3float(int n1, int n2, int n3){
  int i1, i2;
  float ***p;
  p = ( float *** ) calloc( n1, sizeof( float ** ) );
  for ( i1 = 0 ; i1 < n1 ; i1++ ) {
    p[i1] = ( float ** ) calloc( n2, sizeof( float * ) ); 
    for ( i2 = 0 ; i2 < n2 ; i2++ ) {
      p[i1][i2] = ( float * ) calloc( n3, sizeof( float ) );
    }
  }
  return p ;
}

void free3float(int n1, int n2, float ***p){
  int i1, i2;  
  for ( i1 = 0 ; i1 < n1 ; i1++ ) {
    for ( i2 = 0 ; i2 < n2 ; i2++ ) {
      free( p[i1][i2] ) ;
    }
    free( p[i1] ) ;
  }
  free( p ) ;
}

float ****calloc4float(int n1, int n2, int n3, int n4){
  int i1, i2, i3;
  float ****p;
  p = ( float **** ) calloc( n1, sizeof( float *** ) );
  for ( i1 = 0 ; i1 < n1 ; i1++ ) {
    p[i1] = ( float *** ) calloc( n2, sizeof( float ** ) ); 
    for ( i2 = 0 ; i2 < n2 ; i2++ ) {
      p[i1][i2] = ( float ** ) calloc( n3, sizeof( float * ) );
      for ( i3 = 0 ; i3 < n3 ; i3++ ) {
	p[i1][i2][i3] = ( float * ) calloc( n4, sizeof( float ) );
      }	
    }
  }
  return p ;
}

void free4float(int n1, int n2, int n3, float ****p){
  int i1, i2, i3;  
  for ( i1 = 0 ; i1 < n1 ; i1++ ) {
    for ( i2 = 0 ; i2 < n2 ; i2++ ) {
      for ( i3 = 0 ; i3 < n3 ; i3++ ) {
	free( p[i1][i2][i3] ) ;
      }
      free( p[i1][i2] ) ;
    }
    free( p[i1] ) ;
  }
  free( p ) ;
}


int *calloc1int(int n1){
  int *p;
  p = ( int * ) calloc( n1, sizeof( int ) );
  return p ;
}

void free1int(int *p){
  free( p ) ;
}

int **calloc2int(int n1, int n2){
  int i1;
  int **p;
  p = ( int ** ) calloc( n1, sizeof( int * ) );
  for ( i1 = 0 ; i1 < n1 ; i1++ ) {
    p[i1] = ( int * ) calloc( n2, sizeof( int ) ); 
  }
  return p ;
}

void free2int(int n1, int **p){
  int i1 ;  
  for ( i1 = 0 ; i1 < n1 ; i1++ ) {
    free( p[i1] ) ;
  }
  free( p ) ;
}

int ***calloc3int(int n1, int n2, int n3){
  int i1, i2;
  int ***p;
  p = ( int *** ) calloc( n1, sizeof( int ** ) );
  for ( i1 = 0 ; i1 < n1 ; i1++ ) {
    p[i1] = ( int ** ) calloc( n2, sizeof( int * ) ); 
    for ( i2 = 0 ; i2 < n2 ; i2++ ) {
      p[i1][i2] = ( int * ) calloc( n3, sizeof( int ) );
    }
  }
  return p ;
}

void free3int(int n1, int n2, int ***p){
  int i1, i2;  
  for ( i1 = 0 ; i1 < n1 ; i1++ ) {
    for ( i2 = 0 ; i2 < n2 ; i2++ ) {
      free( p[i1][i2] ) ;
    }
    free( p[i1] ) ;
  }
  free( p ) ;
}

int ****calloc4int(int n1, int n2, int n3, int n4){
  int i1, i2, i3;
  int ****p;
  p = ( int **** ) calloc( n1, sizeof( int *** ) );
  for ( i1 = 0 ; i1 < n1 ; i1++ ) {
    p[i1] = ( int *** ) calloc( n2, sizeof( int ** ) ); 
    for ( i2 = 0 ; i2 < n2 ; i2++ ) {
      p[i1][i2] = ( int ** ) calloc( n3, sizeof( int * ) );
      for ( i3 = 0 ; i3 < n3 ; i3++ ) {
	p[i1][i2][i3] = ( int * ) calloc( n4, sizeof( int ) );
      }	
    }
  }
  return p ;
}

void free4int(int n1, int n2, int n3, int ****p){
  int i1, i2, i3;  
  for ( i1 = 0 ; i1 < n1 ; i1++ ) {
    for ( i2 = 0 ; i2 < n2 ; i2++ ) {
      for ( i3 = 0 ; i3 < n3 ; i3++ ) {
	free( p[i1][i2][i3] ) ;
      }
      free( p[i1][i2] ) ;
    }
    free( p[i1] ) ;
  }
  free( p ) ;
}


/*
double complex ****calloc4doublecplx(int n1, int n2, int n3, int n4){
  int i1, i2, i3;
  double complex ****p;
  p = ( double complex **** ) calloc( n1, sizeof( double complex *** ) );
  for ( i1 = 0 ; i1 < n1 ; i1++ ) {
    p[i1] = ( double complex *** ) calloc( n2, sizeof( double complex ** ) ); 
    for ( i2 = 0 ; i2 < n2 ; i2++ ) {
      p[i1][i2] = ( double complex ** ) calloc( n3, sizeof( double complex * ) );
      for ( i3 = 0 ; i3 < n3 ; i3++ ) {
	p[i1][i2][i3] = ( double complex * ) calloc( n4, sizeof( double complex ) );
      }	
    }
  }
  return p ;
}

void free4doublecplx(int n1, int n2, int n3, double complex ****p){
  int i1, i2, i3;  
  for ( i1 = 0 ; i1 < n1 ; i1++ ) {
    for ( i2 = 0 ; i2 < n2 ; i2++ ) {
      for ( i3 = 0 ; i3 < n3 ; i3++ ) {
	free( p[i1][i2][i3] ) ;
      }
      free( p[i1][i2] ) ;
    }
    free( p[i1] ) ;
  }
  free( p ) ;
}
*/