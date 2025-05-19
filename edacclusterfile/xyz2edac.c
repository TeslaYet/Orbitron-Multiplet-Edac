#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#define NMAX 200
void main( int argc, char *argv[] )
{
  char line[NMAX] ;
  int n, i ;
  double x,y,z ;
  FILE *fp ;
  if (argc < 2) { printf("Usage: xyz2edac xyzfile") ; exit(1) ; }
  fp = fopen(argv[1],"r") ; 
  fgets(line,NMAX,fp) ;
  sscanf(line,"%d", &n ) ;
  printf("%4d l(A)\n", n ) ;
  fgets(line,NMAX,fp) ;
  for (i=1;i<=n;i++) {
    fgets(line,NMAX,fp) ;
    printf("%4d %s", i, line ) ;
  }
}
