#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<complex.h>
#define NCHARMAX 4000
#define NFILNAML 120


void main( int argc, char *argv[] )
{
  char inputstring[] = "cluster surface on\nemission angle phi 0 360 181\nemission angle theta 0 80 41\nemission angle window 3\n" ;
  char clusterfile[NFILNAML], outfile[NFILNAML], line[NFILNAML] ;
  char wavestring[NCHARMAX+1] ;  
  int nomega, ngst, nfst, iom, igst, ifst, j, nemitters, *emitternb, *emittertype ;
  double omega, weight, efst ;
  FILE *fp, *fpl, *fpe ;

  printf("edac.clus filename = " ) ; scanf("%s",clusterfile) ;
  printf("nb emitters = " ) ; scanf("%d", &nemitters ) ;
  emitternb   = ( int * ) calloc( nemitters, sizeof( int ) ) ; 
  emittertype = ( int * ) calloc( nemitters, sizeof( int ) ) ;
  printf("Enter emitternb[i] emittertype[i] (i=1..%d): ", nemitters ) ;
  for ( j=0;j<nemitters;j++) scanf("%d%d", &emitternb[j], &emittertype[j] ) ;

  fp  = fopen("rpesalms.edac", "r" ) ;
  fscanf(fp,"%d%d%d%*d",&nomega, &ngst, &nfst ) ; 

  fpl = fopen("logfile","w") ;
  fprintf(fpl,"# %2d %d %3d\n", nomega,  ngst, nfst ) ; 
  
  for ( iom = 0 ; iom < nomega ; iom++ ) { 
    fscanf(fp,"%lf",&omega) ; fprintf(fpl,"%2d %lf\n", iom, omega ) ;
    for ( igst = 0 ; igst < ngst ; igst++ ) { 
      fscanf(fp,"%lf",&weight) ; fprintf(fpl,"   %d %8.6lf\n", igst, weight ) ;
      fgets(line,NCHARMAX,fp )	; 
      for ( ifst = 0 ; ifst < nfst ; ifst++ ) {
	fgets(line,NCHARMAX,fp ) ; sscanf(line,"%lf",&efst) ;
	fprintf(fpl,"   %3d %lf\n", ifst, efst ) ;
	fgets( wavestring,NCHARMAX,fp ) ;
	fpe = fopen("edac.in","w") ;
	fprintf(fpe,"cluster input %s\n", clusterfile ) ;
	sprintf(outfile,"o_%d_%d_%d.ms",iom,igst,ifst) ;
	/* fprintf(fpl,"%2d %d %3d: %lf %lf %lf\n",iom,igst,ifst, omega, weight, efst ) ; */
	fprintf(fpe,"%s", inputstring ) ;
	fprintf(fpe,"emission energy E(eV) %.4lf %.4lf 1\n", omega + efst, omega + efst ) ;
	fprintf(fpe,"emitters %d", nemitters ) ;
	for ( j=0;j<nemitters;j++) fprintf(fpe," %d %d", emitternb[j], emittertype[j] ) ;
	fprintf(fpe,"\nmuffin-tin\norders 1 5\n");
	fprintf(fpe,"V0 E(eV) 10.\nimfp imfp_input\n");
	fprintf(fpe,"wf manual on\nwf manual 2 5 %s", wavestring) ;
	fprintf(fpe,"report emitters\nscan pd %s\n", outfile ) ;
	fprintf(fpe,"report size\nend\n") ;
	fclose(fpe) ;
	system("./run_edac.sh");
      }
    }
  }
  fclose(fpl) ; fclose(fp) ;
}
