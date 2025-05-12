/* typedef's */

/* typedef unsigned int Dualrep ; */
typedef unsigned long long int Dualrep ;


/* define's */

/* for struct Fock */
#define NBITS (8*sizeof(Dualrep)) 
#define ZERO  ((Dualrep)0) 
#define ONE   ((Dualrep)1) 


/* struct's */

struct Fock { Dualrep n, s ; } ;

/* 1-particle operator \sum_{i1,i2} <i2| O1p |i1> c^+_{i2} c_{i1}  */
struct O1p  { int n ; int *i ; struct O1p2 *p ; } ;
struct O1p2 { int n ; int *i ; double *v ; } ;
/* linked list for mat.ele's of O1p. Item members: i1,i2, v == <i2| O1p |i1> */
struct O1plistitem { int i1, i2 ; double v ; struct O1plistitem *next ; } ;

/* 2-particle operator 
   \sum_{i1,i2,i3,i4} <i4,i3| O2p |i1,i2> c^+{i4} c^+_{i3} c_{i2} c_{i1} */
struct O2p  { int n ; int *i ; struct O2p2 *p ; } ;
struct O2p2 { int n ; int *i ; struct O2p3 *p ; } ;
struct O2p3 { int n ; int *i ; struct O2p4 *p ; } ;
struct O2p4 { int n ; int *i ; double *v ; } ;
/* linked list for O2p. Item members: i1,i2,i3,i4, v == <i4,i3| op2p |i1,i2> */
struct O2plistitem { int i1,i2,i3,i4 ; double v ; struct O2plistitem *next ; };

/* one line of a sparse matrix */
struct Spamaline { int n, *j ; double *v ; } ;

/* list of C.I. Slater integrals Rk(ish4,ish3;ish1,ish2) == rmx[k], k=0..nk-1*/
struct CIlistitem { int ish1, ish2, ish3, ish4, nk ; double *rmx ;
                    struct CIlistitem *next ; } ;


/* function's */

/* general */
int sporbindex( int *sorb1sh, int *lsh, int ish, int m, int sigma ) ;
int noverk( int n, int k ) ;
int combinations( int n, int k, int offset, int ***comb ) ;
void readshells( int *pnshells, int **plsh, int **psorb1sh, double **pksish) ;
void readconfs( int nshells, int *lsh, int *pnconfs, int *pnelectrons,
		int ***pocc, int *pnstates ) ;
void makestates( int nshells, int *lsh, int *sorb1sh, int nconfs, 
		 int nelectrons, int **occ, int nstates0, 
		 struct Fock **pstate) ;
int spinorbit2list( int nshells, int *lsh, int *sorb1sh, double *ksish, 
		    struct O1plistitem **ppo1plistitem0 ) ;
int magneticfield2list( int nshells, int *lsh, int *sorb1sh, double hz,
			struct O1plistitem **ppo1plistitem0 ) ;
int exchangefield2list( int nshells, int *lsh, int *sorb1sh, int ish,
			double hz, struct O1plistitem **ppo1plistitem0 ) ;
int exchfieldxz2list( int nshells, int *lsh, int *sorb1sh, int ish, double h, 
		      double theta, struct O1plistitem **ppo1plistitem0 ) ;
int crystalfieldC42list( int nshells, int *lsh, int *sorb1sh, int ish,
                         double ea1, double eb1, double eb2, double ee,
                         struct O1plistitem **ppo1plistitem0 ) ;
int crystalfield2list( int nshells, int *lsh, int *sorb1sh, int ish,
		       double **cfdmat, struct O1plistitem **ppo1plistitem0 ) ;
int coulomb2list( int nshells, int *lsh, int *sorb1sh, 
		  struct CIlistitem *pcilistitem0, 
		  struct O2plistitem **ppo2plistitem0 ) ;
int esh012list( int nshells, int *lsh, int *sorb1sh, double esh0, double esh1,
		struct O1plistitem **ppo1plistitem0 ) ;
int t2gocc2list( int nshells, int *lsh, int *sorb1sh, int ish,
		 struct O1plistitem **ppo1plistitem0 ) ;
int dxyocc2list( int nshells, int *lsh, int *sorb1sh, int ish,
		 struct O1plistitem **ppo1plistitem0 ) ;
int calcham( struct Spamaline **pham, int nstates, struct Fock *state,
	     int nconfs, int nshells, int *lsh, int *sorb1sh, int **occ, 
             double *ksish, double e2p, double e3d, double **cfdmat,
	     double h, double theta ) ;
int solve_dsyev(int n, double* a, double* lambda ) ;
int dipole2list( int nshells, int *lsh, int *sorb1sh, int inish, int q,
                 double **radipmatele, struct O1plistitem **ppitem0 ) ;
int vai2list( int nshells, int *lsh, int *sorb1sh, 
              struct O2plistitem **ppo2plistitem0 ) ;

/* for struct Fock */
void fockcalcs( struct Fock *state ) ;
void fockinitdual( struct Fock *state, Dualrep nin ) ;
void fockinit( struct Fock *state, int norb, int *occorb ) ;
int fockocc( struct Fock *state, int iorb ) ;
int fockcre( struct Fock *state, int iorb ) ;
int fockdes( struct Fock *state, int iorb ) ;
void fockdisplay( struct Fock *state ) ;
void fockdisplays ( struct Fock *state ) ; 
void focklinearsort( int nstates, struct Fock *state ) ;
int findstate( int nstates, struct Fock *state0, struct Fock *statek, int *k );

/* for struct O1plistitem and struct O1p */
void o1plistitemadd( struct O1plistitem **ppitem0, int i1, int i2, double v );
void o1plistdelete( struct O1plistitem *pitem ) ;
struct O1p o1pmake( struct O1plistitem *listp0 ) ;
void o1pprint( struct O1p op1p ) ;
void o1pdelete( struct O1p op1p ) ;
int o1ptospama( int nketbasis, struct Fock *pketbasis, int nbrabasis, 
		struct Fock *pbrabasis, struct O1p op1p, 
		struct Spamaline **ppsml0 ) ;

/* for struct O2plistitem and struct O2p */
void o2plistitemadd( struct O2plistitem **ppitem0,
		int i1, int i2, int i3, int i4, double v ) ;
void o2plistdelete( struct O2plistitem *pitem ) ;
struct O2p o2pmake( struct O2plistitem *listp0 ) ;
void o2pprint( struct O2p op2p ) ;
void o2pdelete( struct O2p op2p ) ;
int o2ptospama( int nketbasis, struct Fock *pketbasis, int nbrabasis, 
		struct Fock *pbrabasis, struct O2p op2p, 
		struct Spamaline **ppsml0 ) ;

/* for struct Spamaline */
int darray2spamaline( struct Spamaline *sml, int n, double *a ) ;
int spamatranspose( struct Spamaline **psm, int nlines, int ncolumns ) ;
void spamalinewrite( struct Spamaline sml ) ;
void spamadelete( struct Spamaline *spama, int nlines ) ;
double spamarealmatele( int nketbasis, double *ketvec, int nbrabasis,
     double *bravec, struct Spamaline *psmline0 ) ;
