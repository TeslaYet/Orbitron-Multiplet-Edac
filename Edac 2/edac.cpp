#ifdef jga_edac            
#else                     
#define jga_edac 1        
                          
#ifdef jga_jga            
#else                     
#define jga_jga 1         
                          
#include <stdio.h>        
#include <stdlib.h>       
#include <string.h>       
#include <time.h>         
#include <math.h>         
                          
#define numero double     
                          
        numero                                        
                                                      
           pi      = 3.141592653589793238462643,      
                                                      
           au_eV   = 27.2113834,         
           a0_au   = 0.529177208,        
           nm      = 10/a0_au,           
           c_au    = 137.03599976;       
                                         
           #define infinite      1.0e20  
           #define infinity      1.0e20  
           #define infinite_int  10000   
           #define infinitesimal 1.0e-20 
int  verbose=0;              
int  relativistic_flag=1;
void relativistic_on (void)  {c_au = 137.03599976;  relativistic_flag=1;}
void relativistic_off(void)  {c_au = 137035.99976;  relativistic_flag=0;}
typedef enum {electrones, fotones, sonido} particle_type_def;
particle_type_def  particle_type=fotones;
FILE *foutput=stdout;        
#define fact_max    100               
#define fact2_max   200               
double  *fact,                        
        *fact2,                       
        *ffact,                       
        *ffact2;                      
int     fact_initialized=0,           
        ffact_initialized=0;          
double factorial(int n)               
{
  double p=1;
  if(n<0) {printf("fact: negative number"); return 1;}
  while(n) p*=n--;
  return p;
}
int init_fact(void)                   
{
  if(fact_initialized)  return 0;
  int n;
  fact   = new double [fact_max+1];       fact_initialized=1;
  fact2  = new double [fact2_max+1];
  fact[0]=fact2[0]=1;
  for(n=1; n<=fact_max; n++)  {fact[n]=n*fact[n-1];  fact2[n]=sqrt(fact[n]);}
  for(n=fact_max+1; n<=fact2_max; n++)  fact2[n]=sqrt(n)*fact2[n-1];
  return 0;
}
int init_ffact(void)                  
{                                     
  if(ffact_initialized)  return 0;
  int n;    init_fact();
  ffact   = new double [fact2_max+1];     ffact_initialized=1;
  ffact2  = new double [fact2_max+1];
  ffact[0]=ffact2[0]=ffact[1]=1;  ffact2[1]=sqrt(pi);
  for(n=2; n<=fact2_max; n++)  {
     ffact[n]=n*ffact[n-2];
     if(n%2==0)  ffact2[n]=fact[n/2-1];
     else        ffact2[n]=(n-2)*(ffact2[n-2]/2);
  }
  return 0;
}
numero combination(int n, int k)      
{
  if(fact_initialized)  return (fact[n]     /(fact[k]     *fact[n-k]     ));
  else                  return (factorial(n)/(factorial(k)*factorial(n-k)));
}
numero sqr(numero a) {return a*a;}                        
int    sqr(int    i) {return i*i;}                        
numero ABS(numero a) {return (a<0)?-a:a;}                 
int    ABS(int    i) {return (i<0)?-i:i;}                 
numero SIGN(numero a) {return (a<0)?-1:((a>0)?1:0);}      
int    SIGN(int    i) {return (i<0)?-1:((i>0)?1:0);}      
int sign_1l(int l)  {if(l%2)  return -1;  return 1;}      
numero *copy(numero *x, int n)
{
  numero *val;  val=new numero [n];
  int i;
  for(i=0; i<n; i++)  val[i]=x[i];
  return val;
}
int *copy(int *x, int n)
{
  int *val;  val=new int [n];
  int i;
  for(i=0; i<n; i++)  val[i]=x[i];
  return val;
}
void cartesian_to_spherical_mu(numero x, numero y, numero z,
                               numero &d, numero &cth, numero &fi)
{
  d=sqrt(x*x+y*y+z*z);                        
  cth=z/d;                                    
  if(!x&&!y) fi=0; else fi=atan2(y,x);        
}
void cartesian_to_spherical(numero x, numero y, numero z,
                            numero &d, numero &th, numero &fi)
{
  cartesian_to_spherical_mu(x,y,z, d,th,fi);  
  if(th>1)  th=0;  else  if(th<-1)  th=pi;    
  th=acos(th);
}
void normalize_angles(numero &th, numero &fi)        
{                                                    
  while(th>pi) th-=2*pi;   while(th<-pi) th+=2*pi;
  if(th<0) {th=-th; fi=pi+fi;}
  while(fi>pi) fi-=2*pi;   while(fi<=-pi) fi+=2*pi;
}
int strcmpC(char *str1, char *str2)
{
  return strcmp(str1,str2);           
}
char read_command[80];                           
void read_name(FILE *f, char *command)
{
  fscanf(f,"%s",command);  char c;  int i;       
                                                 
  if(strlen(command)>1)                          
  if(command[0]=='/' && command[1]=='*') {
    c=1;
    do {
      fscanf(f,"%s",command);  i=strlen(command);
      if(i>1)  if(command[i-2]=='*')  if(command[i-1]=='/')  c=0;
    } while(c);
    read_name(f, command);
  }  else
  if(command[0]=='/' && command[1]=='/') {
    do fscanf(f,"%c",&c); while(c!='\n');
    read_name(f,command);
} }
numero read_numero(char *command)
{
  if(!strcmpC(command,"infinite"))   return infinity;
  if(!strcmpC(command,"infinity"))   return infinity;
  if(!strcmpC(command,"zero"))       return 0;
                                     return atof(command);
}
int read_int(char *command)
{
  if(!strcmpC(command,"infinite"))   return 3*infinite_int;
  if(!strcmpC(command,"infinity"))   return 3*infinite_int;
  if(!strcmpC(command,"all"))        return 2*infinite_int;
  if(!strcmpC(command,"zero"))       return 0;
  if(!strcmpC(command,"off"))        return 0;
  if(!strcmpC(command,"on"))         return 1;
                                     return atoi(command);
}
numero read_numero(FILE *f)
{
  read_name(f,read_command);
  return read_numero(read_command);
}
int read_int(FILE *f)
{
  read_name(f,read_command);
  return read_int(read_command);
}
numero read_numero(void)  {return read_numero(stdout);}
int    read_int(void)     {return read_int(stdout);}
char ftostr_data[50];
char *ftostr(double x, int d)
{
  if(d>0)  if(x-floor(x)==0)  return ftostr(x,0);
  char *val,*v;  val=ftostr_data;  v=val;
  double dec=1e20, eps=1e-12;
  int i=0, j;
  if(x<0) {val[i]='-';  i++;  x=ABS(x);}
  while(floor(x/dec+eps)==0 && dec>1)  dec=dec/10;
  while(dec>=1) {
    j=((int) floor(x/dec+eps));
    x=x-j*dec;  val[i]='0'+j;  i++;  dec=dec/10;
  }
  if(ABS(x)>eps) {
    val[i]='.';  i++;
    while(d) {
      j=int(floor(x/dec+eps));   x=x-j*dec;
      if(d==1 && j<9) if(x/dec>=0.5) j++;
      val[i]='0'+j;  i++;  dec=dec/10;
      d--;
  } }
  val[i]='\0';
  return v;
}
void on_warning(FILE *fout, char *where, char *com1, numero x, char *com2)
  {fprintf(fout,"*** warning in '%s': %s %g %s\n",where,com1,x,com2);}
void on_error(FILE *fout, char *where, char *com1, numero x, char *com2)
  {fprintf(fout,"*** error in '%s': %s %g %s\n",where,com1,x,com2); exit(1);}
void on_warning(FILE *fout, char *where, char *com1, char *com2)
  {fprintf(fout,"*** warning in '%s': %s %s\n",where,com1,com2);}
void on_error(FILE *fout, char *where, char *com1, char *com2)
  {fprintf(fout,"*** error in '%s': %s %s\n",where,com1,com2); exit(1);}
void on_warning(FILE *fout, char *where, char *com)
  {fprintf(fout,"*** warning in '%s': %s\n",where,com);}
void on_error(FILE *fout, char *where, char *com)
  {fprintf(fout,"*** error in '%s': %s\n",where,com); exit(1);}
void on_warning(char *where, char *com1, numero x, char *com2)
  {on_warning(foutput, where,com1,x,com2);}
void on_error(char *where, char *com1, numero x, char *com2)
  {on_error(foutput, where,com1,x,com2);}
void on_warning(char *where, char *com1, char *com2)
  {on_warning(foutput, where,com1,com2);}
void on_error(char *where, char *com1, char *com2)
  {on_error(foutput, where,com1,com2);}
void on_warning(char *where, char *com)  {on_warning(foutput, where,com);}
void on_error(char *where, char *com)  {on_error(foutput, where,com);}
#endif  
#ifdef jga_constants      
#else                     
#define jga_constants 1   
                          
                          
                          
                                         
        numero                           
                                         
           c_mks   = 2.99792458e+8,      
           e_mks   = 1.60217646e-19,     
           hbar_mks= 1.0545716e-34,      
           kb_mks  = 1.3806503e-23,      
           kb_au   = 3.1668151e-06,      
           G_mks   = 6.673e-11,          
           s_mks   = 5.6704e-8,          
           Na      = 6.022142e+23,       
                                         
           au_m    = 0.529177208e-10,    
           au_kg   = 9.1093819e-31,      
           au_s    = 2.4188843e-17,      
           au_v    = 2.18769125e+6,      
                                         
           amu_au  = 1822.88848,         
           mp_au   = 1836.15267,         
           mn_au   = 1838.711873,        
           g_electron=-2.0023193043737,  
                                                      
           ee      = 2.718281828459045235360287,      
           g_Euler = 0.577215664901532860606512;      
                                                      
#endif  
#ifdef jga_files          
#else                     
#define jga_files 1       
                          
                          
FILE *open_file(FILE *fout, char *name, char *mode)
{
  FILE *ff;  ff=NULL;
  if(!strcmpC(name,"inline"))   return stdin;   
  if(!strcmpC(name,"terminal")) return stdout;  
  if(strcmpC(name,"none"))                      
  if((ff=fopen(name,mode))==NULL)
    on_error(fout, "open_file", "cannot open file", name);
  return ff;
}
FILE *open_file(char *name, char *mode)
  {return open_file(foutput,name,mode);}
void close_file(FILE *ff)
  {if(ff!=NULL && ff!=stdout && ff!=stdin) fclose(ff);}
void skip(FILE *fin, int ls)
{
  int i;  char c;
  for(i=0; i<ls; i++)  do fscanf(fin,"%c",&c); while(c!='\n' && !feof(fin));
}
int number_of_columns(int ls, char *name)
{
  FILE *fin;
  int  nc=0;
  char c;
  fin=open_file(name,"r");
  skip(fin,ls);
  if(!feof(fin)) {
               do fscanf(fin,"%c",&c); while(c==' '||c=='\t'||c=='\n');
    do {nc++;  do fscanf(fin,"%c",&c); while(c!=' '&&c!='\t'&&c!='\n');
      if(c!='\n') do fscanf(fin,"%c",&c); while(c==' '||c=='\t');
    } while(c!='\n');
  }
  close_file(fin);
  return nc;
}
int number_of_words(int ls, char *name)
{
  FILE *fin;
  char palabra[80];
  int  nw=-1;       
  fin=open_file(name,"r");
  skip(fin,ls);
  if(!feof(fin))  do {fscanf(fin,"%s",palabra);  nw++;} while(!feof(fin));
  close_file(fin);
  return nw;
}
int number_of_rows(int ls, char *name)
{
  int nc=number_of_columns(ls,name);
  int nw=number_of_words(ls,name);
  if(nw%nc!=0)  on_warning("number_of_rows", "unstructured file");
  return  nw/nc;
}
int number_of_chars(int ls, char *name)
{
  FILE *fin;
  int  nc=0;
  char c;
  fin=open_file(name,"r");
  skip(fin,ls);
  if(!feof(fin))  do {fscanf(fin,"%c",&c);  nc++;} while(!feof(fin));
  close_file(fin);
  return nc;
}
int number_of_columns(char *name) {return number_of_columns(0,name);}
int number_of_words  (char *name) {return number_of_words  (0,name);}
int number_of_rows   (char *name) {return number_of_rows   (0,name);}
int number_of_chars  (char *name) {return number_of_chars  (0,name);}
void limit_values(int ls, char *name, int &nc,
                  numero *vmin, numero *vmax, numero *vave)
{
  FILE *fin;
  int  nr=number_of_rows(ls,name);  nc=number_of_columns(ls,name);
  int  i,j;
  float  v;
  for(j=0; j<nc; j++)  {vmin[j]=infinity; vmax[j]=-infinity; vave[j]=0;}
  fin=open_file(name,"r");
  skip(fin,ls);
  for(i=0; i<nr; i++)
  for(j=0; j<nc; j++) {
     fscanf(fin,"%f",&v);  vave[j]+=v;
     if(v>vmax[j])  vmax[j]=v;
     if(v<vmin[j])  vmin[j]=v; 
  }
  close_file(fin);
  for(j=0; j<nc; j++)  vave[j]=vave[j]/nr;
}
#endif  
#ifdef jga_time           
#else                     
#define jga_time 1        
                          
                          
#define time_seconds  0.000001  
numero  time_initial=clock(),   
        time_partial=0,         
        time_total=0;           
void time_update(void)
{
  numero time_temp=clock();
  time_partial+=time_temp-time_initial;
  time_initial=time_temp;
}
void time_print(FILE *fout)
{
  time_update();  time_total+=time_partial;
  fprintf(fout, "--- Time:  total=%4.2f sec, partial= %4.2f sec.\n",
          time_total*time_seconds,
          time_partial*time_seconds);
  time_partial=0;
}
void time_print(void) {time_print(foutput);}
#endif  
#ifdef jga_cg             
#else                     
#define jga_cg 1          
                          
                          
double Clebsch_Gordan(int a, int b, int c, int alpha, int beta, int gamma)
{
  double  delabc,term1,term2=0,term3;   int  nu=0;
  init_fact();
  if(alpha + beta != gamma || a<0 || b<0 || c<0 ||
     ABS(alpha)>a || ABS(beta)>b || ABS(gamma)>c ||
     c<ABS(a-b) || a+b<c)                           return 0;
  delabc = ((fact2[a+b-c]/fact2[a+b+c+1]) * fact2[a+c-b]) * fact2[b+c-a];
  term1 = sqrt((2*c+1))*fact2[a+alpha]*fact2[a-alpha]*fact2[b+beta];
  term1 *= fact2[b-beta]*fact2[c+gamma]*fact2[c-gamma];
  while((nu <= (a-alpha)) && (nu <= (b+beta)) && (nu <= (a+b-c))) {
    if((c-b+alpha+nu >= 0) && (c-a-beta+nu >= 0)) {
      term3 = ((((((term1/fact[a-alpha-nu])/fact[c-b+alpha+nu])
                         /fact[b+beta -nu])/fact[c-a-beta +nu])
                         /fact[nu])/fact[a+b-c-nu]);
      term2 += (((nu<0?-nu:nu)&1)?-term3:term3);
    }
    nu++;
  }
  return delabc*term2;
}
double Clebsch_Gordan(int a, int b, int c, int alpha, int beta)
  {return  Clebsch_Gordan(a,b,c,alpha,beta,alpha+beta);}
double Gaunt(int a, int b, int c, int alpha, int beta)
{
  int gamma=alpha+beta;
  if(a<0 || b<0 || c<0 ||
     ABS(alpha)>a || ABS(beta)>b || ABS(gamma)>c ||
     c<ABS(a-b) || a+b<c)                           return 0;
  return  sqrt((2*a+1)*(2*b+1)/(4*pi*(2*c+1)))
            *Clebsch_Gordan(a,b,c,0,0,0)
            *Clebsch_Gordan(a,b,c,alpha,beta,gamma);
}
double Gaunt(int a, int b, int c, int alpha, int beta, int gamma)
{
  if(alpha + beta != gamma)   return 0;
  return Gaunt(a,b,c,alpha,beta);
}
double CGspin1_2(int j1, int j3, int m2, int m3)
{
  if(j1<0 || j3<0 || m3>j3 || m3<-j3 || m3-m2>j1 || m3-m2<-j1)  return 0;
  if(j3==j1+1)  if(m2==1)    return   sqrt((j1+m3+1.0)/(2*j1+2.0));  else
                if(m2==-1)   return   sqrt((j1-m3+1.0)/(2*j1+2.0));
  if(j3==j1-1)  if(m2==1)    return  -sqrt((j1-m3+1.0)/(2*j1+2.0));  else
                if(m2==-1)   return   sqrt((j1+m3+1.0)/(2*j1+2.0));
                             return   0;
}
double CGspin1(int j1, int j3, int m2, int m3)
{
  if(j3==j1+2)  if(m2==2)  return   sqrt(  (j1+m3)*(j1+m3+2.0)
                                         /((2*j1+2.0)*(2*j1+4.0)));  else
                if(m2==0)  return   sqrt(  (j1-m3+2.0)*(j1+m3+2.0)
                                         /((2*j1+2.0)*(j1+2.0)));    else
                if(m2==-2) return   sqrt(  (j1-m3)*(j1-m3+2.0)
                                         /((2*j1+2.0)*(2*j1+4.0)));
  if(j3==j1)    if(m2==2)  return  -sqrt(  (j1+m3)*(j1-m3+2.0)
                                         /( 2*j1*(j1+2.0)));         else
                if(m2==0)  return   m3/sqrt(j1*(j1+2.0));            else
                if(m2==-2) return   sqrt(  (j1-m3)*(j1+m3+2.0)
                                         /( 2*j1*(j1+2.0)));
  if(j3==j1-2)  if(m2==2)  return   sqrt(  (j1-m3)*(j1-m3+2.0)
                                         /( 2*j1*(2*j1+2.0)));       else
                if(m2==0)  return  -sqrt(  (j1-m3)*(j1+m3)
                                         /( j1*(2*j1+2.0)));         else
                if(m2==-2) return   sqrt(  (j1+m3+2.0)*(j1+m3)
                                         /( 2*j1*(2*j1+2.0)));
                           return   0;
}
int k_lj(int l, int j)  {if(j==2*l+1)  return -l-1;  return l;}
void lj_k(int &l, int &j, int k)
{
  if(k<0)  {l=-k-1; j=2*l+1;}  else  {l=k; j=2*l-1;}
}
double Gaunt_spin(int k1, int k2, int mj1, int mj2, int i)
{
  int j1,l1,j2,l2;
  lj_k(l1,j1,k1);  lj_k(l2,j2,k2);  if(l1!=l2) return 0;  l1=2*l1;
  if(ABS(mj1)>j1 || ABS(mj2)>j2)  return 0;
  if(i==-1)  if(mj1!=mj2-2)  return 0;
             else            return      CGspin1_2(l1,j1,-1,mj1)
                                       * CGspin1_2(l1,j2, 1,mj2) / sqrt(2);
  if(i== 0)  if(mj1!=mj2  )  return 0;
             else            return (    CGspin1_2(l1,j1, 1,mj1)
                                       * CGspin1_2(l1,j2, 1,mj1)
                                     -   CGspin1_2(l1,j1,-1,mj1)
                                       * CGspin1_2(l1,j2,-1,mj1) )/2;
  if(i== 1)  if(mj1!=mj2+2)  return 0;
             else            return      CGspin1_2(l1,j1, 1,mj1)
                                       * CGspin1_2(l1,j2,-1,mj2) / sqrt(2);
  return 0;  
}
#endif  
#ifdef jga_complex        
#else                     
#define jga_complex 1     
                          
                          
class complex {
  friend complex operator+(numero   x, const complex &v);
  friend complex operator+(const complex &u, numero   x);
  friend complex operator+(const complex &u, const complex &v);
  friend complex operator-(const complex &v);
  friend complex operator-(numero   x, const complex &v);
  friend complex operator-(const complex &u, numero   x);
  friend complex operator-(const complex &u, const complex &v);
  friend complex operator*(numero   x, const complex &v);
  friend complex operator*(const complex &u, numero   x);
  friend complex operator*(const complex &u, const complex &v);
  friend complex operator/(numero   x, const complex &v);
  friend complex operator/(const complex &u, numero   x);
  friend complex operator/(const complex &u, const complex &v);
  friend numero  real(const complex &v);
  friend numero  imag(const complex &v);
  friend numero  mod(const complex  &v);
  friend numero  mod2(const complex &v);
  friend complex conj(const complex &v);
  friend complex exp(const complex  &v);
  friend complex pow(const complex  &v, int m);
  friend complex pow(const complex  &a, numero x);
  friend complex sqrt(const complex &v);
  friend complex sqr(const complex &v);
  friend complex log(const complex  &v);
  friend complex llog(const complex  &v);
  friend complex sin(const complex  &v);
  friend complex cos(const complex  &v);
  friend complex atan(const complex  &v);
  friend complex sinh(const complex &v);
  friend complex cosh(const complex &v);
  friend complex euler(numero a, numero b);
  friend numero  arg(complex a);
public:
  numero &real(void);
  numero &imag(void);
  complex(void) {}
  complex(numero re);
  complex(int re);
  complex(numero re, numero im);
  complex(const complex &v);
  complex &operator=(numero    x);
  complex &operator=(int       x);
  complex &operator=(const complex  &v);
  complex &operator+=(numero   x);
  complex &operator+=(const complex &v);
  complex &operator-=(numero   x);
  complex &operator-=(const complex &v);
private:
  numero re, im;
};
inline numero &complex::real(void)
       {return re;}
inline numero &complex::imag(void)
       {return im;}
inline complex::complex(numero r)
       {re=r; im=0.0;}
inline complex::complex(int    r)
       {re=r; im=0.0;}
inline complex::complex(numero r, numero i)
       {re=r; im=i;}
inline complex::complex(const complex &v)
       {re=v.re; im=v.im;}
inline complex &complex::operator=(numero    x)
       {re=x; im=0.0; return *this;}
inline complex &complex::operator=(int       x)
       {re=x; im=0.0; return *this;}
inline complex &complex::operator=(const complex  &v)
       {re=v.re; im=v.im; return *this;}
inline complex &complex::operator+=(numero   x)
       {re+=x; return *this;}
inline complex &complex::operator+=(const complex &v)
       {re+=v.re; im+=v.im; return *this;}
inline complex &complex::operator-=(numero   x)
       {re-=x; return *this;}
inline complex &complex::operator-=(const complex &v)
       {re-=v.re; im-=v.im; return *this;}
inline complex operator+(numero   x, const complex &v)
       {return complex(x+v.re,v.im);}
inline complex operator+(const complex &u, numero   x)
       {return complex(u.re+x,u.im);}
inline complex operator+(const complex &u, const complex &v)
       {return complex(u.re+v.re,u.im+v.im);}
inline complex operator-(const complex &v)
       {return complex(-v.re,-v.im);}
inline complex operator-(numero   x, const complex &v)
       {return complex(x-v.re,-v.im);}
inline complex operator-(const complex &u, numero   x)
       {return complex(u.re-x,u.im);}
inline complex operator-(const complex &u, const complex &v)
       {return complex(u.re-v.re,u.im-v.im);}
inline complex operator*(numero   x, const complex &v)
       {return complex(x*v.re,x*v.im);}
inline complex operator*(const complex &u, numero   x)
       {return complex(u.re*x,u.im*x);}
inline complex operator*(const complex &u, const complex &v)
       {return complex(u.re*v.re-u.im*v.im,u.re*v.im+u.im*v.re);}
inline complex operator/(numero   x, const complex &v)
       {numero m=v.re*v.re+v.im*v.im;
        return complex(x*v.re/m,-x*v.im/m);}
inline complex operator/(const complex &u, numero   x)
       {return complex(u.re/x,u.im/x);}
inline complex operator/(const complex &u, const complex &v)
       {numero m=v.re*v.re+v.im*v.im;
        return complex((u.re*v.re+u.im*v.im)/m,(u.im*v.re-u.re*v.im)/m);}
inline numero real(const complex &v)
       {return v.re;}
inline numero imag(const complex &v)
       {return v.im;}
inline numero mod(const complex  &v)
       {return sqrt(v.re*v.re+v.im*v.im);}
inline numero mod2(const complex &v)
       {return v.re*v.re+v.im*v.im;}
inline complex conj(const complex &v)
       {return complex(v.re,-v.im);}
inline complex exp(const complex  &v)
       {return exp(v.re)*complex(cos(v.im),sin(v.im));}
inline complex sqr(const complex &v)
       {return complex(v.re*v.re-v.im*v.im, 2*(v.re*v.im));}
inline complex euler(numero a, numero b)
       {return complex(a*cos(b), a*sin(b));}
inline complex sin(const complex &v)
       {return complex(sin(v.re)*cosh(v.im), cos(v.re)*sinh(v.im));}
inline complex cos(const complex &v)
       {return complex(cos(v.re)*cosh(v.im),-sin(v.re)*sinh(v.im));}
inline complex sinh(const complex &v)
       {return (exp(v)-exp(-v))/2;}
inline complex cosh(const complex &v)
       {return (exp(v)+exp(-v))/2;}
complex pow(const complex &v, int m)
{
  if(m<0)   return pow(1/v,-m);
  if(m==0)  return complex(1,0);
  if(m==1)  return v;
  if(m==2)  return v*v;
  if(m%2)   return v*sqr(pow(v,m/2));
            return   sqr(pow(v,m/2));
}
complex pow(const complex  &a, numero x)
{
  if(x<0)   return pow(1/a,-x);
  if(x==0)  return complex(1,0);
            return euler(pow(mod(a),x),x*arg(a));
}
complex sqrt(const complex &v)  
{
  if(!v.im)  if(v.re<0)  return complex(0,sqrt(-v.re));
             else        return complex(sqrt(v.re), 0);
  numero a;
  if(v.re>0) {
    a=sqrt((v.re+sqrt(v.re*v.re+v.im*v.im))/2);
    return complex(a,v.im/(2*a));
  }  else {
    a=sqrt((-v.re+sqrt(v.re*v.re+v.im*v.im))/2);
    if(v.im>=0)  return complex(v.im/(2*a),a);
    else         return complex(-v.im/(2*a),-a);
} }
numero arg(complex a)
{
  if(!a.re && !a.im)  return 0;
  return atan2(a.im,a.re);
}
complex log(const complex  &v)
{
  complex b=0;
  if(!v.re && !v.im)  printf(" Error: log(0)\n");
  numero va=mod(v);
  if(va<1e-20)  return llog(v);
  b.re=log(va);
  if(!v.re)  if(v.im>0) b.im=pi/2; else b.im=-pi/2;
  else {
    b.im=atan(fabs(v.im/v.re));
    if(v.im<0 && v.re>0)  b.im=-b.im;        else
    if(v.im<0 && v.re<0)  b.im= b.im +  pi;  else
    if(v.im>0 && v.re<0)  b.im=-b.im +  pi;
    while(b.im>pi)   b.im-=2*pi;
    while(b.im<-pi)  b.im+=2*pi;
  }
  return b;
}
complex llog(const complex  &v)   
{
  numero magnitude=fabs(v.re)+fabs(v.im);
  complex b=0, w=v/magnitude;
  if(!v.re && !v.im)  printf(" Error: llog(0)\n");
  b.re=log(mod(w))+log(magnitude);
  if(!v.re)  if(v.im>0) b.im=pi/2; else b.im=-pi/2;
  else {
    b.im=atan(fabs(v.im/v.re));
    if(v.im<0 && v.re>0)  b.im=-b.im;        else
    if(v.im<0 && v.re<0)  b.im= b.im +  pi;  else
    if(v.im>0 && v.re<0)  b.im=-b.im +  pi;
    while(b.im>pi)   b.im-=2*pi;
    while(b.im<-pi)  b.im+=2*pi;
  }
  return b;
}
complex atan(const complex &v)
{
  exit(0);
  return complex(cos(v.re)*cosh(v.im),-sin(v.re)*sinh(v.im));
}
complex i_c(0,1);
complex i_l(int l)
{
  l=l%4;  if(l<0) l+=4;
  if(l==0)  return complex( 1, 0);
  if(l==1)  return complex( 0, 1);
  if(l==2)  return complex(-1, 0);
  if(l==3)  return complex( 0,-1);
  return 0;
}
complex *copy(complex *x, int n)
{
  complex *val;  val=new complex [n];
  int i;
  for(i=0; i<n; i++)  val[i]=x[i];
  return val;
}
#endif  
#ifdef jga_qfint          
#else                     
#define jga_qfint 1       
                          
                          
numero *qfint_xv, *qfint_wv;
int qfint_mmax=0;
int qfint_init(int mmax)
{
  if(mmax<=qfint_mmax) return 0;
  numero  pi2=2*pi, theta;
  int  i,j,m;
  if(qfint_mmax>0) {delete [] qfint_xv;  delete [] qfint_wv;}
  qfint_xv=new numero [mmax+2];
  qfint_wv=new numero [mmax+2];
  j=1;  m=4;  qfint_mmax=mmax;
  do {
    for(i=1; i<m; i++)  if(i%2 || m==4) {
      theta = i * (pi2 / m);
      qfint_xv [j]= i/((numero) m) - sin (theta) / pi2;
      qfint_wv [j]=  1  - cos (theta);
      j++;
    }
    m=2*m;
  } while (m<=mmax || m<=8);
  return 0;
}
numero qfint(numero a, numero b, numero func(numero x), numero err, int mmax)
{
  numero  esm,sigmam;   int  i,j,m;
  numero   eta, e;   
  if(mmax>qfint_mmax)  qfint_init(mmax);
  if(a==b)  return 0;
  m=4;  j=1;  esm=0;
  do {
    sigmam=0;
    for(i=1; i<m; i++)  if(i%2 || m==4) {
      eta=(b-a)*qfint_xv[j]+a;
      sigmam+=qfint_wv[j]*func(eta);
      j++;
    }
    sigmam=2*(b-a)/m*sigmam;
    e     =(esm-sigmam)/2;   if(e<0) e=-e;
    esm   =0.5*(esm+sigmam);
    m=2*m;
  } while ((e>err && m<=mmax)  || m<=8);
  return  esm;
}
#ifdef jga_complex
complex qfint(numero a, numero b, complex func(numero x), numero err, int mmax)
{
  complex  esm,sigmam;   int  i,j,m;
  numero   eta, e;   
  if(mmax>qfint_mmax)  qfint_init(mmax);
  if(a==b)  return 0;
  m=4;  j=1;  esm=0;
  do {
    sigmam=0;
    for(i=1; i<m; i++)  if(i%2 || m==4) {
      eta=(b-a)*qfint_xv[j]+a;
      sigmam+=qfint_wv[j]*func(eta);
      j++;
    }
    sigmam=2*(b-a)/m*sigmam;
    e     =mod(esm-sigmam)/2;
    esm   =0.5*(esm+sigmam);
    m=2*m;
  } while ((e>err && m<=mmax)  || m<=8);
  return  esm;
}
#endif  
#endif  
#ifndef jga_bessel        
#define jga_bessel 1      
                          
                          
numero chebev(numero a, numero b, numero c[], int m, numero x)
{
  numero d=0,dd=0,sv,y,y2;  int j;
  if((x-a)*(x-b)>0)  on_warning("chebev", "x not in range");
  y=(2*x-a-b)/(b-a);  y2=2*y;
  for(j=m; j>=1; j--) {sv=d;  d=y2*d-dd+c[j];  dd=sv;}
  return y*d-dd+c[0]/2;
}
void beschb (numero x, numero &g1, numero &g2, numero &gamp, numero &gamm)
{
  numero xx,c1[7],c2[8];
  c1[0] = -1.142022680371172;  c1[1] = 0.006516511267076;
  c1[2] = 0.000308709017308;   c1[3] = -3.470626964e-6;
  c1[4] = 6.943764e-9;         c1[5] = 3.678e-11;
  c1[6] = -1.36e-13;
  c2[0] = 1.843740587300906;   c2[1] = -0.076852840844786;
  c2[2] = 0.001271927136655;   c2[3] = -4.971736704e-6;
  c2[4] = -3.3126120e-8;       c2[5] = 2.42310e-10;
  c2[6] = -1.7e-13;            c2[7] = -1.0e-15;
  xx=8*x*x-1;
  g1=chebev(-1,1,c1,6,xx);     g2=chebev(-1,1,c2,7,xx);
  gamp=g2-x*g1;                gamm=g2+x*g1;
}
int besselJY(numero xnu, numero x, numero &rj,  numero &ry,
             numero &rjp, numero &ryp)
{
  int i,isign=1,l,nl, maxit=100000, salir;
  numero xmin=2.0, eps=1.0e-10, fpmin=1.0e-30,
	 a,b,br,bi,c,cr,ci,d=0,del,del1,den,di,dlr,dli,
       	 dr,e,f,fct,fct2,fct3,ff,gam,gam1,gam2,gammi,gampl,h,
      	 p,pimu,pimu2,q,r,rjl,rjl1,rjmu,rjp1,rjpl,rjtemp,ry1,
       	 rymu,rymup,rytemp,sum,sum1,temp,w,x2,xi,xi2,xmu,xmu2;
  if(x<=0)  on_error("besselJY", "wrong argument x =", x, "");
  if(xnu<0)  on_error("besselJY", "wrong order n =", xnu, "");
  if(x<xmin)    nl=int(xnu+0.5);
  else {nl=int(xnu-x+1.5); nl=0>nl?0:nl;}
  xmu=xnu-nl;  xmu2=xmu*xmu;  xi=1/x;  xi2=2*xi;  w=xi2/pi;  h=xnu*xi;
  if(h<fpmin) h=fpmin;     c=h;  b=xi2*xnu;
  for(i=1, salir=1; i<=maxit && salir; i++) {
    b+=xi2;  d=b-d;             if(fabs(d)<fpmin) d=fpmin;
    c=b-1/c;                    if(fabs(c)<fpmin) c=fpmin;
    d=1/d;  del=c*d;  h=del*h;  if(d<0) isign=-isign;
    if(fabs(del-1)<eps) salir=0; else salir=1;
  }
  if(salir)  on_warning("besselJY", "x too large");
  rjl1=rjl=isign*fpmin;   rjp1=rjpl=h*rjl;  fct=xnu*xi;
  for(l=nl; l; l--) {rjtemp=fct*rjl+rjpl;
                     fct-=xi;  rjpl=fct*rjtemp-rjl;  rjl=rjtemp;}
  if(rjl==0.0) rjl=eps;
  f=rjpl/rjl;
  if(x<xmin) {
    x2=x/2;  pimu=pi*xmu;
    if(fabs(pimu)<eps)  fct=1;  else fct=pimu/sin(pimu);
    d=-log(x2);  e=xmu*d;
    if(fabs(e)<eps)     fct2=1; else fct2=sinh(e)/e;
    beschb(xmu,gam1,gam2,gampl,gammi);
    ff=2/pi*fct*(gam1*cosh(e)+gam2*fct2*d);
    e=exp(e);   p=e/(gampl*pi);
    q=1/(e*pi*gammi);
    pimu2=pimu/2;
    if(fabs(pimu2)<eps) fct3=1; else fct3=sin(pimu2)/pimu2;
    r=pi*pimu2*fct3*fct3;
    c=1;  d=-x2*x2;  sum=ff+r*q;  sum1=p;
    for(i=1, salir=1; i<=maxit && salir; i++) {
      ff=(i*ff+p+q)/(i*i-xmu2);   c*=d/i;
      p/=(i-xmu);   q/=(i+xmu);   del=c*(ff+r*q);
      sum+=del;   del1=c*p-i*del;
      sum1+=del1;
      if(fabs(del)<(1+fabs(sum))*eps) salir=0; else salir=1;
    }
    if(salir)  on_warning("besselJY", "series failed to converge");
    rymu=-sum;  ry1=-sum1*xi2;  rymup=xmu*xi*rymu-ry1;
    rjmu=w/(rymup-f*rymu);
  } else {
    a=0.25-xmu2;  p=-xi/2;  q=1;  br=2*x;  bi=2;
    fct=a*xi/(p*p+q*q);  cr=br+q*fct;  ci=bi+p*fct;
    den=br*br+bi*bi;  dr=br/den;  di=-bi/den;  dlr=cr*dr-ci*di;
    dli=cr*di+ci*dr;  temp=p*dlr-q*dli;  q=p*dli+q*dlr;
    p=temp;
    for(i=2, salir=1; i<=maxit && salir; i++) {
      a+=2*(i-1);  bi+=2;  dr=a*dr+br;  di=a*di+bi;
      if(fabs(dr)+fabs(di)<fpmin) dr=fpmin;
      fct=a/(cr*cr+ci*ci);  cr=br+cr*fct;  ci=bi-ci*fct;
      if(fabs(cr)+fabs(ci)<fpmin) cr=fpmin;
      den=dr*dr+di*di;  dr=dr/den;  di=-di/den;
      dlr=cr*dr-ci*di;  dli=cr*di+ci*dr;  temp=p*dlr-q*dli;
      q=p*dli+q*dlr;  p=temp;
      if(fabs(dlr-1)+fabs(dli)<eps) salir=0; else salir=1;
    }
    if(salir)  on_warning("besselJY", "cf2 failed");
    gam=(p-f)/q;
    rjmu=sqrt(w/((p-f)*gam+q));
    if(rjl<0) rjmu=-rjmu;
    rymu=rjmu*gam;
    rymup=rymu*(p+q/gam);
    ry1=xmu*xi*rymu-rymup;
  }
  fct=rjmu/rjl;
  rj=rjl1*fct;
  rjp=rjp1*fct;
  for(i=1; i<=nl; i++) {
    rytemp=(xmu+i)*xi2*ry1-rymu;
    rymu=ry1;  ry1=rytemp;
  }
  ry=rymu;  ryp=xnu*xi*rymu-ry1;
  return 0;
}
void besselIK(numero xnu, numero x, numero &ri, numero &rk, numero &rip,
              numero &rkp)
{
  int maxit=100000;
  int i,l,nl;
  numero xmin=2.0, eps=1.0e-10, fpmin=1.0e-30,
         a,a1,b,c,d,del,del1,delh,dels,e,f,fact,fact2,ff,gam1,gam2,
         gammi,gampl,h,p,pimu,q,q1,q2,qnew,ril,ril1,rimu,rip1,ripl,
         ritemp,rk1,rkmu,rkmup,rktemp,s,sum,sum1,x2,xi,xi2,xmu,xmu2;
  
  if(x<=0)  on_error("besselIK", "wrong argument x =", x, "");
  if(xnu<0)  on_error("besselIK", "wrong order n =", xnu, "");
  nl=(int)(xnu+0.5);
  xmu=xnu-nl;
  xmu2=xmu*xmu;
  xi=1.0/x;
  xi2=2.0*xi;
  h=xnu*xi;
  if (h < fpmin) h=fpmin;
  b=xi2*xnu;
  d=0.0;
  c=h;
  for (i=1;i<=maxit;i++) {
    b += xi2;
    d=1.0/(b+d);
    c=b+1.0/c;
    del=c*d;
    h=del*h;
    if (fabs(del-1.0) < eps) break;
  }
  if (i > maxit)  on_warning("besselIK", "x too large: try asymptotic expan.");
  ril=fpmin;
  ripl=h*ril;
  ril1=ril;
  rip1=ripl;
  fact=xnu*xi;
  for (l=nl;l>=1;l--) {
    ritemp=fact*ril+ripl;
    fact -= xi;
    ripl=fact*ritemp+ril;
    ril=ritemp;
  }
  f=ripl/ril;
  if (x < xmin) {
    x2=0.5*x;
    pimu=pi*xmu;
    fact = (fabs(pimu) < eps ? 1.0 : pimu/sin(pimu));
    d = -log(x2);
    e=xmu*d;
    fact2 = (fabs(e) < eps ? 1.0 : sinh(e)/e);
    beschb(xmu,gam1,gam2,gampl,gammi);
    ff=fact*(gam1*cosh(e)+gam2*fact2*d);
    sum=ff;
    e=exp(e);
    p=0.5*e/gampl;
    q=0.5/(e*gammi);
    c=1.0;
    d=x2*x2;
    sum1=p;
    for (i=1;i<=maxit;i++) {
      ff=(i*ff+p+q)/(i*i-xmu2);
      c *= (d/i);
      p /= (i-xmu);
      q /= (i+xmu);
      del=c*ff;
      sum += del;
      del1=c*(p-i*ff);
      sum1 += del1;
      if (fabs(del) < fabs(sum)*eps) break;
    }
    if (i > maxit) on_warning("besselIK", "series failed to converge");
    rkmu=sum;
    rk1=sum1*xi2;
  } else {
    b=2.0*(1.0+x);
    d=1.0/b;
    h=delh=d;
    q1=0.0;
    q2=1.0;
    a1=0.25-xmu2;
    q=c=a1;
    a = -a1;
    s=1.0+q*delh;
    for (i=2;i<=maxit;i++) {
      a -= 2*(i-1);
      c = -a*c/i;
      qnew=(q1-b*q2)/a;
      q1=q2;
      q2=qnew;
      q += c*qnew;
      b += 2.0;
      d=1.0/(b+a*d);
      delh=(b*d-1.0)*delh;
      h += delh;
      dels=q*delh;
      s += dels;
      if (fabs(dels/s) < eps) break;
    }
    if (i > maxit)  on_warning("besselIK", "cf2 failed");
    h=a1*h;
    rkmu=sqrt(pi/(2.0*x))*exp(-x)/s;
    rk1=rkmu*(xmu+x+0.5-h)*xi;
  }
  rkmup=xmu*xi*rkmu-rk1;
  rimu=xi/(f*rkmu-rkmup);
  ri=(rimu*ril1)/ril;
  rip=(rimu*rip1)/ril;
  for (i=1;i<=nl;i++) {
    rktemp=(xmu+i)*xi2*rk1+rkmu;
    rkmu=rk1;
    rk1=rktemp;
  }
  rk=rkmu;
  rkp=xnu*xi*rkmu-rk1;
}
numero besselJ0(numero x)
{
  numero  ax,xx,z,y,ans1,ans2;
  if(x==0)  return 1;
  if(x<8 && x>-8)  {
      y = sqr(x);
      ans1 = 57568490574.0+y*(-13362590354.0+y*(651619640.7
         +y*(-11214424.18+y*(77392.33017+y*(-184.9052456)))));
      ans2 = 57568490411.0+y*(1029532985.0+y*(9494680.718
         +y*(59272.64853+y*(267.8532712+y*1.0))));
      return ans1/ans2;
  }
  ax=x<0?-x:x; z = 8.0/ax; y = sqr(z); xx = ax-0.785398164;
  ans1 = 1.0+y*(-0.1098628627e-2+y*(0.2734510407e-4
     +y*(-0.2073370639e-5+y*0.2093887211e-6)));
  ans2 = -0.1562499995e-1+y*(0.1430488765e-3
     +y*(-0.6911147651e-5+y*(0.7621095161e-6
     -y*0.934945152e-7)));
  return  sqrt(0.636619772/ax)*(cos(xx)*ans1-z*sin(xx)*ans2);
}
numero besselJ1(numero x)
{
  numero  ax,xx,z,y,ans1,ans2;
  if (x==0)  return 0;
  if (x<8 && x>-8)  {
    y = sqr(x);
    ans1 = x*(72362614232.0+y*(-7895059235.0+y*(242396853.1
       +y*(-2972611.439+y*(15704.48260+y*(-30.16036606))))));
    ans2 = 144725228442.0+y*(2300535178.0+y*(18583304.74
       +y*(99447.43394+y*(376.9991397+y*1.0))));
    return  ans1/ans2;
  }
  ax = x<0?-x:x; z = 8.0/ax; y = sqr(z); xx = ax-2.356194491;
  ans1 = 1.0+y*(0.183105e-2+y*(-0.3516396496e-4
     +y*(0.2457520174e-5+y*(-0.240337019e-6))));
  ans2 = 0.04687499995+y*(-0.2002690873e-3
     +y*(0.8449199096e-5+y*(-0.88228987e-6+y*0.105787412e-6)));
  return  sqrt(0.636619772/ax)*(cos(xx)*ans1
          -z*sin(xx)*ans2)*(x>0?1:-1);
}
numero besselJJ(int n, numero x)
{
  if(n<0) if(n%2) return -besselJJ(-n,x); else return besselJJ(-n,x);
  if(ABS(x)<1e-12) if(n==0) return 1; else return 0;
  numero J,Y,Jp,Yp;  besselJY(n,x,J,Y,Jp,Yp);
  return J;
}
numero besselYY(int n, numero x)
{
  if(n<0) if(n%2) return -besselJJ(-n,x); else return besselJJ(-n,x);
  if(ABS(x)<1e-12) if(n==0) return 1; else return 0;
  numero J,Y,Jp,Yp;  besselJY(n,x,J,Y,Jp,Yp);
  return Y;
}
numero besselJp(int n, numero x)
{
  if(n<0) if(n%2) return -besselJp(-n,x); else return besselJp(-n,x);
  if(ABS(x)<1e-12) if(n==1)  return 0.5;  else
                   if(n==-1) return -0.5; else return 0;
  numero J,Y,Jp,Yp;  besselJY(n,x,J,Y,Jp,Yp);
  return Jp;
}
numero besselJ(int n, numero x)
{
  int iacc=40, j,jsum,m;
  numero  bigno=1.0e10, bigni=1.0e-10,
          bj,bjm,bjp,sum,tox,ans=0;
  if(ABS(x)<1e-12) if(n==0) return 1; else return 0;
  if(n==0)  return  besselJ0(x);
  if(n==1)  return  besselJ1(x);
  if(n<0)   if(n%2)  return -besselJ(-n,x);  else  return besselJ(-n,x);
  tox = 2.0/x;
  if(x>n)  {
    bjm=besselJ0(x);
    bj=besselJ1(x);
    for(j=1; j<n; j++) {
      bjp=j*tox*bj-bjm;
      bjm=bj;
      bj=bjp;
    }
    return  bj;
  }
  m= 2*((n+ int(floor(sqrt(1.0*(iacc*n)))) ) / 2);
  jsum=0;  sum=bjp=0;  bj=1;
  for(j=m; j>=1; j--) {
    bjm=j*tox*bj-bjp;
    bjp=bj;
    bj=bjm;
    if(bj>bigno || bj<-bigno) {
      bj=bj*bigni;
      bjp=bjp*bigni;
      ans=ans*bigni;
      sum=sum*bigni;
    }
    if(jsum)  sum = sum+bj;
    jsum = 1-jsum;
    if(j==n)  ans=bjp;
  }
  sum = 2.0*sum-bj;
  return  ans/sum;
}
void besseljy(int l, numero x, numero &j, numero &y)
{
  numero jp,yp;
  numero factor=sqrt(pi/(2*x));
  besselJY(l+0.5,x,j,y,jp,yp);
  j*=factor;  y*=factor;
}
numero besselj(int l, numero x)
{
  if(ABS(x)<1e-10)  if(l==0) return 1.0;     else
                    if(l==1) return x/3;     else
                    if(l==2) return x*x/15;  else  return 0;
  numero jl,yl;
  besseljy(l,x,jl,yl);
  return jl;
}
numero bessely(int l, numero x)
{
  numero jl,yl;
  besseljy(l,x,jl,yl);
  return yl;
}
void besseljyp(int l, numero x, numero &jp, numero &yp)
{
  numero j,y;
  numero factor=sqrt(pi/(2*x));
  besselJY(l+0.5,x,j,y,jp,yp);
  jp=factor*(jp-j/(2*x));  yp=factor*(yp-y/(2*x));
}
numero besselI0(numero x)
{
  numero  ax,y;
  if(-3.75<x && x<3.75) {
    y=x/3.75;  y=y*y;
    return  1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492+y*
                  (0.2659732+y*(0.360768e-1+y*0.45813e-2)))));
  } else {
    ax=(x<0)?-x:x;  y=3.75/ax;
    return  (exp(ax)/sqrt(ax))*(0.39894228+y*(0.1328592e-1
              +y*(0.225319e-2+y*(-0.157565e-2+y*(0.916281e-2
              +y*(-0.2057706e-1+y*(0.2635537e-1+y*(-0.1647633e-1
              +y*0.392377e-2))))))));
  }
}
numero besselI1(numero x)
{
  numero  ax,y,ans;
  if(-3.75<x && x<3.75) {
    y=x/3.75;  y=y*y;
    return  x*(0.5+y*(0.87890594+y*(0.51498869+y*(0.15084934
               +y*(0.2658733e-1+y*(0.301532e-2+y*0.32411e-3))))));
  } else {
    ax=(x<0)?-x:x;  y=3.75/ax;
    ans=0.2282967e-1+y*(-0.2895312e-1+y*(0.1787654e-1-y*0.420059e-2));
    ans=0.39894228+y*(-0.3988024e-1+y*(-0.362018e-2
                  +y*(0.163801e-2+y*(-0.1031555e-1+y*ans))));
    return  (exp(ax)/sqrt(ax))*ans;
  }
}
numero besselK0(numero x)
{
  numero  y;
  if(x<=2) {
    y=x*x/4;
    return  -log(x/2.0)*besselI0(x)+(-0.57721566+y*(0.42278420
               +y*(0.23069756+y*(0.3488590e-1+y*(0.262698e-2
               +y*(0.10750e-3+y*0.74e-5))))));
  } else {
    y=2/x;
    return  (exp(-x)/sqrt(x))*(1.25331414+y*(-0.7832358e-1
               +y*(0.2189568e-1+y*(-0.1062446e-1+y*(0.587872e-2
               +y*(-0.251540e-2+y*0.53208e-3))))));
  }
}
numero besselK1(numero x)
{
  numero  y;
  if(x<=2) {
      y=x*x/4;
      return  log(x/2.0)*besselI1(x)+(1.0/x)*(1.0+y*(0.15443144
                  +y*(-0.67278579+y*(-0.18156897+y*(-0.1919402e-1
                  +y*(-0.110404e-2+y*(-0.4686e-4)))))));
  } else {
      y=2/x;
      return  (exp(-x)/sqrt(x))*(1.25331414+y*(0.23498619
                +y*(-0.3655620e-1+y*(0.1504268e-1+y*(-0.780353e-2
                 +y*(0.325614e-2+y*(-0.68245e-3)))))));
  }
}
numero besselI(int n, numero x)
{
  if(ABS(x)<1e-20)  if(n==0)  return 1.0;  else  return 0.0;
  int     iacc=40, j,m;
  numero  bigno=1.0e10,
          bigni=1.0e-10,
          bi=1, bim, bip=0, tox=2/x, ans=0;
  if(n==0)  return  besselI0(x);
  if(n==1)  return  besselI1(x);
  if(n<0)   n=-n;
  m=2*(n+int(floor(sqrt(iacc*n))));
  for(j=m; j>=1; j--) {
    bim=bip+j*tox*bi;    bip=bi;    bi=bim;
    if(bi<-bigno || bigno<bi) {
         ans=ans*bigni;
         bi=bi*bigni;
         bip=bip*bigni;
    }
    if(j==n)  ans=bip;
  }
  return  ans*besselI0(x)/bi;
}
numero besselK(int n, numero x)
{
  numero  tox,bkp,bkm,bk;
  int  j;
  if(n==0)  return  besselK0(x);
  if(n==1)  return  besselK1(x);
  if(n<0)   n=-n;
  tox=2/x;
  bkm=besselK0(x);
  bk=besselK1(x);
  for(j=1; j<n; j++) {
    bkp=bkm+j*tox*bk;    bkm=bk;    bk=bkp;
  }
  return  bk;
}
numero besselIp(int n, numero x)
{
  if(n==0)  return besselI(1,x);
  return  (besselI(n-1,x)+besselI(n+1,x))/2;
}
numero besselKp(int n, numero x)
{
  if(n==0)  return -besselK(1,x);
  return -(besselK(n-1,x)+besselK(n+1,x))/2;
}
#ifdef jga_complex
int besselIKexp_flag=1;
complex besselI0(complex x)
{
  complex  y,eax,deax;  if(real(x)<0) on_error("besselI0", "real(x)<0");
  if(besselIKexp_flag) {eax=exp(x); deax=1;} else {eax=1; deax=exp(-x);}
  if(mod(x)<3.75) {
    y=x/3.75;  y=y*y;
    return  (1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492
                +y*(0.2659732+y*(0.360768e-1+y*0.45813e-2))))))*deax;
  } else {
    y=3.75/x;
    return  (eax/sqrt(x))*(0.39894228+y*(0.1328592e-1
              +y*(0.225319e-2+y*(-0.157565e-2+y*(0.916281e-2
              +y*(-0.2057706e-1+y*(0.2635537e-1+y*(-0.1647633e-1
              +y*0.392377e-2))))))));
} }
complex besselI1(complex x)
{
  complex  y,ans,eax,deax;  if(real(x)<0) on_error("besselI1", "real(x)<0");
  if(besselIKexp_flag) {eax=exp(x); deax=1;} else {eax=1; deax=exp(-x);}
  if(mod(x)<3.75) {
    y=x/3.75;  y=y*y;
    return  (x*(0.5+y*(0.87890594+y*(0.51498869+y*(0.15084934
                   +y*(0.2658733e-1+y*(0.301532e-2+y*0.32411e-3)))))))*deax;
  } else {
    y=3.75/x;
    ans=0.2282967e-1+y*(-0.2895312e-1+y*(0.1787654e-1-y*0.420059e-2));
    ans=0.39894228+y*(-0.3988024e-1+y*(-0.362018e-2
                  +y*(0.163801e-2+y*(-0.1031555e-1+y*ans))));
    return  (eax/sqrt(x))*ans;
} }
complex besselK0(complex x)
{
  complex  y,eax,deax;
  if(besselIKexp_flag) {eax=exp(-x); deax=1;} else {eax=1; deax=exp(x);}
  if(mod(x)<=2) {
    y=x*x/4;
    return  (-log(x/2.0)*besselI0(x)*deax+(-0.57721566+y*(0.42278420
               +y*(0.23069756+y*(0.3488590e-1+y*(0.262698e-2
               +y*(0.10750e-3+y*0.74e-5)))))))*deax;
  } else {
    y=2/x;
    return  eax/sqrt(x)*(1.25331414+y*(-0.7832358e-1
               +y*(0.2189568e-1+y*(-0.1062446e-1+y*(0.587872e-2
               +y*(-0.251540e-2+y*0.53208e-3))))));
} }
complex besselK1(complex x)
{
  complex  y,eax,deax;
  if(besselIKexp_flag) {eax=exp(-x); deax=1;} else {eax=1; deax=exp(x);}
  if(mod(x)<=2) {
      y=x*x/4;
      return  (log(x/2.0)*besselI1(x)*deax+(1.0/x)*(1.0+y*(0.15443144
                  +y*(-0.67278579+y*(-0.18156897+y*(-0.1919402e-1
                  +y*(-0.110404e-2+y*(-0.4686e-4))))))))*deax;
  } else {
    y=2/x;
    return  eax/sqrt(x)*(1.25331414+y*(0.23498619
              +y*(-0.3655620e-1+y*(0.1504268e-1+y*(-0.780353e-2
              +y*(0.325614e-2+y*(-0.68245e-3)))))));
} }
complex besselI(int n, complex x)
{
  if(mod(x)<1e-20)  if(n==0)  return 1.0;  else  return 0.0;
  int     iacc=40, j,m;
  numero   bigno=1.0e10, bigni=1.0e-10;
  complex   bi(1,0), bim, bip(0,0), tox, ans=0;   tox=2/x;
  if(n<0)   n=-n;
  if(n==0)  return  besselI0(x);
  if(n==1)  return  besselI1(x);
  m=2*(n+int(floor(sqrt(iacc*n))));
  for(j=m; j>=1; j--) {
    bim=bip+j*tox*bi;    bip=bi;    bi=bim;
    if(bigno<mod(bi)) {
         ans=ans*bigni;
         bi=bi*bigni;
         bip=bip*bigni;
    }
    if(j==n)  ans=bip;
  }
  return  ans*besselI0(x)/bi;
}
complex besselK(int n, complex x)
{
  complex  tox,bkp,bkm,bk;
  int  j;
  if(n<0)   n=-n;
  if(n==0)  return  besselK0(x);
  if(n==1)  return  besselK1(x);
  tox=2/x;
  bkm=besselK0(x);
  bk=besselK1(x);
  for(j=1; j<n; j++) {
    bkp=bkm+j*tox*bk;    bkm=bk;    bk=bkp;
  }
  return  bk;
}
complex besselIp(int n, complex x)
{
  if(n==0)  return besselI(1,x);
  return  (besselI(n-1,x)+besselI(n+1,x))/2;
}
complex besselKp(int n, complex x)
{
  if(n==0)  return -besselK(1,x);
  return -(besselK(n-1,x)+besselK(n+1,x))/2;
}
complex besselh(int l, numero x)
{
  numero jl,yl;
  besseljy(l,x,jl,yl);
  return i_c*jl-yl;
}
complex besselhp(int l, numero x)
{
  numero j,y,jp,yp;
  numero factor=sqrt(pi/(2*x));
  besselJY(l+0.5,x,j,y,jp,yp);
  return factor*((i_c*jp-yp)-(i_c*j-y)/(2*x));
}
complex bessely(int l, complex z)
{
  if(l==0)  return -cos(z)/z;
  if(l==1)  return -(cos(z)/z+sin(z))/z;
  
  int n;
  complex yl_2, yl_1=bessely(0,z), yl=bessely(1,z);    
  for(n=2; n<=l; n++) {yl_2=yl_1;  yl_1=yl;  yl=(n+n-1)/z*yl_1-yl_2;}
  return yl;
}
complex besselyp(int l, complex z)
{
  if(l==0) return -bessely(1,z);
  else     return  bessely(l-1,z)-((l+1)/z)*bessely(l,z);
}
complex besselilh(int l, complex z)
{                                        
  complex val=1, fct=1, fct1=i_c/(2*z);  
  int s;
  init_fact();
  for(s=1; s<=l; s++) {
    fct=fct*fct1/s;
    val=val+fct*fact[l+s]/fact[l-s];
  }
  return  exp(i_c*z)/z * val;
}
complex besselilhp(int l, complex z)
{                                        
  if(l>0) return i_c*besselilh(l-1,z)-(l+1)/z*besselilh(l,z);
          return i_c*besselilh(l+1,z)+    l/z*besselilh(l,z);
}
#ifdef jga_qfint
complex besselj_int_z;
int     besselj_int_l;
complex besselj_int(numero th)
{
  return cos(besselj_int_z*cos(th)) * pow(sin(th),2*besselj_int_l+1.0);
}
complex besselj(int l, complex z)
{
  if(l==0)  if(!real(z) && !imag(z))  return complex(1,0);
            else                      return sin(z)/z;
  if(l==1)  if(!real(z) && !imag(z))  return complex(0,0);
            else                      return (sin(z)/z-cos(z))/z;
  complex coef(1,0);   int i;
  besselj_int_z=z;
  besselj_int_l=l;
  for(i=1; i<=l; i++)  coef=coef*z/(2*i);
  return coef * qfint(0,pi/2,besselj_int,0.00001,2048);
}
complex besseljp(int l, complex z)
{
  if(l==0) return -besselj(1,z);
  else     return  besselj(l-1,z)-((l+1)/z)*besselj(l,z);
}
complex besselh(int l, complex z)
{
  if(ABS(imag(z))<infinitesimal)  return besselh(l, real(z));
  return i_c*besselj(l,z) -bessely(l,z);
}
complex besselhp (int l, complex z)
{
  if(ABS(imag(z))<infinitesimal)  return besselhp(l, real(z));
  return i_c*besseljp(l,z)-besselyp(l,z);
}
complex besselh1 (int l, complex z)  {return besselj(l,z) +i_c*bessely(l,z);}
complex besselh1p(int l, complex z)  {return besseljp(l,z)+i_c*besselyp(l,z);}
complex besselh2 (int l, complex z)  {return besselj(l,z) -i_c*bessely(l,z);}
complex besselh2p(int l, complex z)  {return besseljp(l,z)-i_c*besselyp(l,z);}
void besseljh(int l, numero x, numero &j, numero &jp, complex &h, complex &hp)
{
  numero y,yp, x2=2*x;
  numero factor=sqrt(pi/x2);
  besselJY(l+0.5,x,j,y,jp,yp);
  j=factor*j;         y=factor*y;         h=i_c*j-y;      
  jp=factor*jp-j/x2;  yp=factor*yp-y/x2;  hp=i_c*jp-yp;
}
void besseljh(int l, complex z, complex &j, complex &jp,
                                complex &h, complex &hp)
{
  complex y,yp, hh,hhp, b=complex(0,imag(z));
  numero  rj,rjp,rjj,rjjp, x=real(z);
  if(ABS(imag(z))<1e-2*ABS(real(z)) && ABS(imag(z))<1e-2) {
    besseljh(l,   real(z), rj, rjp, h, hp);
    besseljh(l+1, real(z), rjj,rjjp,hh,hhp);
    j=rj+b*rjp;  h=h+b*hp;
    jp=rjp+b*(l/x*(rjp-rj/x)-rjjp);
    hp= hp+b*(l/x*( hp- h/x)- hhp);
  } else {
    on_warning("besseljh", "complex z");
    j=besselj(l,z);  jp=besseljp(l,z);
    y=bessely(l,z);  yp=besselyp(l,z);
    h=i_c*j-y;       hp=i_c*jp-yp;
} }
#endif  
#endif  
#endif  
#ifdef jga_spline         
#else                     
#define jga_spline 1      
                          
                          
class spline {
public:
  numero *x, *y, *y2, a,b;         
  int n, interpol;
  spline(void) {n=0;}
  void alloc(int na, int type_interpol);
  void alloc(int na) {alloc(na,0);}
  void put(int na, numero xa, numero ya) {x[na]=xa;  y[na]=ya;}
  void add(int na, numero xa, numero ya) {x[na]+=xa;  y[na]+=ya;}
  void prod(numero v) {int i;  for(i=0; i<n; i++) y[i]=y[i]*v;}
  int  init(numero yp0, numero ypn_1);
  void init(void) {if(n>1) init(infinity,infinity);}
  void init(numero *xx, numero *yy, int nn);
  void init(numero *xx, numero *yy, int nn, int type_interpol);
  numero integ(numero aa, numero bb, int l);
  numero integ(numero aa, numero bb) {return integ(aa,bb,0);}
  numero val(numero xx);
  numero der(numero xx);
  void free(void);
  
  void find(numero xx, int &ii, int &jj);  
  numero integ(int k, numero aa, numero bb, int l);
  numero dpow(numero aa, numero bb, int l);
};
void spline::free(void)
{
  if(n>0)  {delete [] x;  delete [] y;}
  if(n>=4 && interpol==0)  delete [] y2;
  n=0;
}
void spline::alloc(int na, int type_interpol)
{
  free();  n=na;  interpol=type_interpol;
  if(n>0)  {
    x=new numero [n];  y=new numero [n];
    for(na=0; na<n; na++) {x[na]=0; y[na]=0;}
  }
  if(n>=4 && interpol==0) y2=new numero [n];
}
int spline::init(numero yp0, numero ypn_1)
{
    int  i,k;   numero  p,qn,sig,un;   numero  *u;
    if(n==0)  return 0;                  
    if(n==1)  {y[0]=yp0;  return 0;}     
    a=x[0];  b=x[n-1];
    if(interpol)  return 0;
    u=new numero [n];
    if(yp0>0.99*infinity)  y2[0]=u[0]=0;
    else  { y2[0]=-0.5;
            u[0]=3/(x[1]-x[0])*((y[1]-y[0])/(x[1]-x[0])-yp0);}
    for(i=1; i<n-1; i++) {
      sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
      p=sig*y2[i-1]+2;
      y2[i]=(sig-1)/p;
      u[i]=(y[i+1]-y[i])/(x[i+1]-x[i])
            -(y[i]-y[i-1])/(x[i]-x[i-1]);
      u[i]=(6*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
    }
    if(ypn_1>0.99*infinity)  qn=un=0;
    else  {qn=0.5;
           un=3.0/(x[n-1]-x[n-2])*(ypn_1-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));}
    y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1);
    for(k=n-2; k>=0; k--)  y2[k]=y2[k]*y2[k+1]+u[k];
    delete [] u;   return 0;
}
numero spline::val(numero xx)
{
    if(n==0)  return 0;
    if(n==1)  return y[0];
    if(n<4 && interpol==0)  on_error("spline::val", "less than 4 points");
    int  klo,khi;   numero  h,bb,aa;
    find(xx, klo,khi);
    h=x[khi]-x[klo];
    if(h==0)  on_error("spline::val", "vanishing spacing");
    if(interpol)  return y[klo] + (y[khi]-y[klo])*(xx-x[klo])/h;
    aa=(x[khi]-xx)/h;
    bb=(xx-x[klo])/h;
    return aa*y[klo]+bb*y[khi]+((aa*aa*aa-aa)*y2[klo]
                               +(bb*bb*bb-bb)*y2[khi])*h*h/6;
}
numero spline::der(numero xx)
{
    if(n<=1)  return 0;
    if(n<4 && interpol==0)  on_error("spline::der", "less than 4 points");
    int  klo,khi;   numero  h,bb,aa,ap,bp;
    find(xx, klo,khi);
    h=x[khi]-x[klo];
    if(h==0)  on_error("spline::der", "vanishing spacing");
    if(interpol)  return (y[khi]-y[klo])/h;
    aa=(x[khi]-xx)/h;    ap=-1/h;
    bb=(xx-x[klo])/h;    bp= 1/h;
    return ap*y[klo]+bp*y[khi]+
           ((3*aa*aa*ap-ap)*y2[klo]+(3*bb*bb*bp-bp)*y2[khi])*h*h/6;
}
void spline::init(numero *xx, numero *yy, int nn, int type_interpol)
{
  free();  n=nn;  interpol=type_interpol;
  if(type_interpol==0)  y2=new numero [n];
  x=copy(xx,n);  y=copy(yy,n);
  init();
}
void spline::init(numero *xx, numero *yy, int nn) {init(xx,yy,nn,0);}
void spline::find(numero xx, int &klo, int &khi)
{
  int s,k;
  if(x[0]>x[n-1])  s=0;  else  s=1;
  klo=0;   khi=n-1;
  while(khi-klo>1) {
    k=(khi+klo)/2;
    if((x[k]>xx&&s)||(x[k]<xx&&s==0))  khi=k;  else  klo=k;
} }
numero spline::dpow(numero aa, numero bb, int l)
{
  if(l!=-1)  return (pow(bb,l+1.0)-pow(aa,l+1.0))/(l+1);
  else       return log(bb/aa);
}
numero spline::integ(int k, numero aa, numero bb, int l)
{
  numero x0=x[k], x1=x[k+1];
  numero h=x1-x0;
  numero I0=dpow(aa,bb,l);
  numero I1=dpow(aa,bb,l+1);
  if(interpol)  return y[k]*I0+(y[k+1]-y[k])*(I1-x[k]*I0)/h;
  numero I2=dpow(aa,bb,l+2);
  numero I3=dpow(aa,bb,l+3);
  numero h3=h*h*h;
  numero A1=(x1*I0-I1)/h;
  numero B1=(I1-x0*I0)/h;
  numero A3=((( I0*x1-3*I1)*x1+3*I2)*x1-I3)/h3;
  numero B3=(((-I0*x0+3*I1)*x0-3*I2)*x0+I3)/h3;
  return A1*y[k]+B1*y[k+1]+((A3-A1)*y2[k]+(B3-B1)*y2[k+1])*h*h/6;
}
numero spline::integ(numero aa, numero bb, int l)
{
  if(n==0) return 0;
  if(n==1) return y[0]*dpow(aa,bb,l);
  if(n<4 && interpol==0) on_error("spline::integ", "less than 4 points");
  int i, kalo,kahi,kblo,kbhi;  find(aa,kalo,kahi);  find(bb,kblo,kbhi);
  numero  val=integ(kalo,aa,x[kahi],l);
  for(i=kahi; i<kblo; i++)  val+=integ(i,x[i],x[i+1],l);
  return  val+integ(kblo,x[kblo],bb,l);
}
#ifdef jga_complex
class splinec {
public:
  numero *x, a,b;         
  complex *y, *y2;
  int n, interpol;
  splinec(void) {n=0;}
  void alloc(int na, int type_interpol);
  void alloc(int na) {alloc(na,0);}
  void put(int na, numero xa, complex ya) {x[na]=xa;  y[na]=ya;}
  void add(int na, numero xa, complex ya) {x[na]+=xa;  y[na]+=ya;}
  void prod(complex v) {int i;  for(i=0; i<n; i++) y[i]=y[i]*v;}
  int init(complex yp0, complex ypn_1);
  void init(void) {if(n>1) init(infinity,infinity);}
  void init(numero *xx, complex *yy, int nn);
  complex val(numero xx);
  complex der(numero xx);
  void free(void);
};
void splinec::free(void)
{
  if(n>0)  {delete [] x;  delete [] y;}
  if(n>=4 && interpol==0)  delete [] y2;
  n=0;
}
void splinec::alloc(int na, int type_interpol)
{
  free();  n=na;  interpol=type_interpol;
  if(n>0)  {
    x= new numero [n];  y= new complex [n];
    for(na=0; na<n; na++) {x[na]=0; y[na]=0;}
  }
  if(n>=4 && interpol==0) y2=new complex [n];
}
int splinec::init(complex yp0, complex ypn_1)
{
    int  i,k;   complex  p,qn,sig,un;   complex  *u;
    if(n==0)  return 0;                  
    if(n==1)  {y[0]=yp0;  return 0;}     
    a=x[0];  b=x[n-1];
    if(interpol)  return 0;
    u=new complex [n];
    if(real(yp0)>0.99*infinity)  y2[0]=u[0]=0;
    else  { y2[0]=-0.5;
            u[0]=3/(x[1]-x[0])*((y[1]-y[0])/(x[1]-x[0])-yp0);}
    for(i=1; i<n-1; i++) {
      sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
      p=sig*y2[i-1]+2;
      y2[i]=(sig-1)/p;
      u[i]=(y[i+1]-y[i])/(x[i+1]-x[i])
            -(y[i]-y[i-1])/(x[i]-x[i-1]);
      u[i]=(6*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
    }
    if(real(ypn_1)>0.99*infinity)  qn=un=0;
    else  {qn=0.5;
           un=3.0/(x[n-1]-x[n-2])*(ypn_1-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));}
    y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1);
    for(k=n-2; k>=0; k--)  y2[k]=y2[k]*y2[k+1]+u[k];
    delete [] u;   return 0;
}
complex splinec::val(numero xx)
{
    if(n==0)  return 0;
    if(n==1)  return y[0];
    if(n<4 && interpol==0)  on_error("splinec::val", "less than 4 points");
    int  klo,khi,k,s;   numero  h,bb,aa;
    if(x[0]>x[n-1])  s=0;  else  s=1;
    klo=0;   khi=n-1;
    while(khi-klo>1) {
      k=(khi+klo)/2;
      if((x[k]>xx&&s)||(x[k]<xx&&s==0))  khi=k;  else  klo=k;
    }
    h=x[khi]-x[klo];
    if(h==0)  on_error("splinec::val", "vanishing spacing");
    if(interpol)  return y[klo] + (y[khi]-y[klo])*(xx-x[klo])/h;
    aa=(x[khi]-xx)/h;
    bb=(xx-x[klo])/h;
    return aa*y[klo]+bb*y[khi]+((aa*aa*aa-aa)*y2[klo]
                               +(bb*bb*bb-bb)*y2[khi])*h*h/6;
}
complex splinec::der(numero xx)
{
    if(n<=1)  return 0;
    if(n<4 && interpol==0)  on_error("splinec::der", "less than 4 points");
    int  klo,khi,k,s;   numero  h,bb,aa,ap,bp;
    if(x[0]>x[n-1])  s=0;  else  s=1;
    klo=0;   khi=n-1;
    while(khi-klo>1) {
      k=(khi+klo)/2;
      if((x[k]>xx&&s)||(x[k]<xx&&s==0))  khi=k;  else  klo=k;
    }
    h=x[khi]-x[klo];
    if(h==0)  on_error("splinec::der", "vanishing spacing");
    if(interpol)  return (y[khi]-y[klo])/h;
    aa=(x[khi]-xx)/h;   ap=-1/h;
    bb=(xx-x[klo])/h;   bp= 1/h;
    return ap*y[klo]+bp*y[khi]+
           ((3*aa*aa*ap-ap)*y2[klo]+(3*bb*bb*bp-bp)*y2[khi])*h*h/6;
}
void splinec::init(numero *xx, complex *yy, int nn)
{
  free();  n=nn;  interpol=0;
  x=copy(xx,n);  y=copy(yy,n);  y2=new complex [n];
  init();
}
#endif  
#endif  
#ifdef jga_legendre       
#else                     
#define jga_legendre 1    
                          
                          
numero legendre(int l, int m, numero x)
{
  if(m<0 || m>l || x<-1 || 1<x) {
    printf(" Error: Bad argument in legendre: l=%1d  m=%1d  x=%f\n",l,m,x);
    return 0;
  }
  int  i,ll;
  numero  pmm=1, somx2, fct=1, pmmp1, pll;
  if(m>0) {
    somx2=sqrt((1-x)*(1+x));
    for(i=1; i<=m; i++, fct+=2)  pmm=pmm*fct*somx2;
  }
  if(l==m)    return pmm;
  pmmp1=x*(2*m+1)*pmm;
  if(l>m+1)  for(ll=m+2; ll<=l; ll++) {
               pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m);
               pmm=pmmp1;  pmmp1=pll;
             }
  return  pmmp1;
}
#ifdef jga_complex
complex legendre(int l, int m, complex x)
{
  if(m<0 || m>l) {
    printf(" Error: Bad argument in legendre: l=%1d  m=%1d\n",l,m);
    return complex(0,0);
  }
  int  i,ll;
  complex  pmm=1, somx2, fct=1, pmmp1, pll;
  if(m>0) {
    somx2=sqrt((1-x)*(1+x));
    for(i=1; i<=m; i++, fct+=2)  pmm=pmm*fct*somx2;
  }
  if(l==m)    return pmm;
  pmmp1=x*(2*m+1)*pmm;
  if(l>m+1)  for(ll=m+2; ll<=l; ll++) {
               pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m);
               pmm=pmmp1;  pmmp1=pll;
             }
  return  pmmp1;
}
#endif  
#endif  
#ifdef jga_Ylm            
#else                     
#define jga_Ylm 1         
                          
                          
numero Ylma(int l, int m)
{
  numero a;  int mm=ABS(m);
  a=sqrt((l+l+1)*fact[l-mm]/(4*pi*fact[l+mm]));
  if(m>0)  if(m%2)  return -a;
                    return  a;
}
numero Ylmmu0(int l, int m, numero cth)
  {return  legendre(l,ABS(m),cth) * Ylma(l,m);}
numero Ylm0(int l, int m, numero th)
{
  while(th> pi) th-=2*pi;
  while(th<-pi) th+=2*pi;
  return Ylmmu0(l,m,cos(th));
}
complex Ylmmu(int l, int m, numero cth, numero fi)
  {return Ylmmu0(l,m,cth) * euler(1,m*fi);}
complex Ylm(int l, int m, numero th, numero fi)
  {return Ylm0(l,m,th) * euler(1,m*fi);}
complex Ylmmu(int l, int m, complex cth, numero fi)
  {return  legendre(l,ABS(m),cth) * Ylma(l,m) * euler(1,m*fi);}
#endif  
#ifdef jga_matrix         
#else                     
#define jga_matrix 1      
                          
                          
class matrix {
public:
  int n,m,nm, sel;                        
  complex *p;                             
  matrix(void) {n=m=nm=sel=0;}
  void free(void)  {if(nm) {delete [] p;}  n=m=nm=0;}
  matrix(int i, int j)  {nm=0;  alloc(i,j);}
  matrix(int i) {nm=0;  alloc(i,1);}
  matrix(const matrix &A);
  ~matrix(void) {free();}
  void alloc(int i, int j);
  void alloc(int i)  {alloc(i,1);}
  void alloc_diagonal(int i) {alloc(i,1);  sel=1;}
  void a(int i, int j, complex x) {p[i*m+j]=x;}
  void a(int i, int j, numero x) {p[i*m+j]=x;}
  void a(int i, complex x) {p[i]=x;}
  void a(int i, numero x) {p[i]=x;}
  void add(int i, int j, complex x) {p[i*m+j]+=x;}
  void add(int i, int j, numero x) {p[i*m+j]+=x;}
  void add(int i, complex x) {p[i]+=x;}
  void add(int i, numero x) {p[i]+=x;}
  void set_column(int j, matrix &B);
  void set_row(int i, matrix &B);
  complex a(int i, int j);
  complex a(int i) {if(i>-1 && i<nm)  return p[i];  else  return complex(0,0);}
  matrix &operator=(complex x);
  matrix &operator=(numero x);
  matrix &operator=(const matrix &A);
  matrix &operator++(void);                
  matrix &operator--(void);                
  complex inv(void);                       
  complex gaussj(const matrix &Mb, int n, int m, int opt);
  complex gaussj(const matrix &Mb, int n, int m)   
    {return gaussj(Mb,n,m,0);}
  void print(numero val);
  void print(void) {print(0);}
  friend matrix transpose(matrix A);
  friend matrix conj(matrix A);
  friend matrix operator+(complex x, matrix A);
  friend matrix operator+(numero x, matrix A) {return complex(x,0)+A;};
  friend matrix operator+(matrix A, complex x) {return x+A;}
  friend matrix operator+(matrix A, numero x) {return x+A;}
  friend matrix operator+(matrix A, const matrix &B);
  friend matrix operator%(matrix A, matrix  B);
  friend matrix operator-(matrix A);
  friend matrix operator-(complex x, matrix A);
  friend matrix operator-(numero x, matrix A) {return complex(x,0)-A;};
  friend matrix operator-(matrix A, complex x) {return (-x)+A;}
  friend matrix operator-(matrix A, numero x) {return (-x)+A;}
  friend matrix operator-(matrix A, const matrix &B);
  friend matrix operator*(complex x, matrix A);
  friend matrix operator*(numero x, matrix A);
  friend matrix operator*(matrix A, complex x) {return x*A;}
  friend matrix operator*(matrix A, numero x) {return x*A;}
  friend matrix operator*(const matrix &A, const matrix &B);
  friend matrix operator*(numero *r, matrix A);
  friend matrix operator*(const matrix &A, numero *r) {return r*A;}
  friend matrix operator%(numero *r, matrix A);
  friend matrix operator%(matrix A, numero *r);
  friend complex prod(matrix A, matrix B);
  friend matrix operator/(complex x, matrix A) {return x*(--A);}
  friend matrix operator/(numero x, matrix A) {return x*(--A);}
  friend matrix operator/(matrix A, complex x) {return (1/x)*A;}
  friend matrix operator/(matrix A, numero x) {return (1/x)*A;}
  friend matrix operator/(matrix A, matrix B) {return A*(--B);}
  friend complex trace(const matrix &A);
  friend complex sum(const matrix &A);
  friend numero mod2(const matrix &A);
  friend numero mod(const matrix &A) {return sqrt(mod2(A));}
  friend int dim(const matrix &A) {return A.nm;}
  friend complex det(const matrix &A);
  friend complex lndet(const matrix &A);
  friend int dimn(const matrix &A) {return A.n;}
  friend int dimm(const matrix &A) {return A.m;}
  friend matrix neg_off_diag(matrix A);
};
matrix::matrix(const matrix &A)
{
  int ij;
  nm=0;  if(A.nm) alloc(A.n,A.m);  sel=A.sel;
  for(ij=0; ij<nm; ij++)  p[ij]=A.p[ij];
}
complex matrix::a(int i, int j)
{
  if(i>-1 && i<n) {
    if(sel==1)  if(i==j)          return p[i];
    if(sel==0)  if (j>-1 && j<m)  return p[i*m+j];
  }
  return 0;
}
void matrix::alloc(int i, int j)
{
  if(nm)  delete [] p;   n=i; m=j; nm=i*j;  sel=0;
  if(nm) {p=new complex [nm];  for(i=0; i<nm; i++) p[i]=0;}
}
void matrix::set_column(int j, matrix &B)
{
  if((B.n!=1 && B.m!=1) || B.nm!=n)
    on_error("matrix::set_column", "incompatible dimensions");
  int i;
  complex *q,*r;  q=p+j;  r=B.p;
  for(i=0; i<n; i++, q+=m, r++) *q=*r;
}
void matrix::set_row(int i, matrix &B)
{
  if((B.n!=1 && B.m!=1) || B.nm!=m)
    on_error("matrix::set_row", "incompatible dimensions");
  int j;
  complex *q,*r;  q=p+i*m;  r=B.p;
  for(j=0; j<m; j++, q++, r++) *q=*r;
}
matrix operator+(complex x, matrix A)
{
  int i,j,ij;
  if(A.sel==0) {
    if(A.n!=A.m)  on_error("matrix","A is not square in x+A");
    for(i=ij=0; i<A.n; i++)  for(j=0; j<A.m; j++,ij++)  if(i==j)  A.p[ij]+=x;
  } else
  if(A.sel==1)  for(i=0; i<A.n; i++)  A.p[i]+=x;
  return A;
}
matrix operator+(matrix A, const matrix &B)
{
  int i,j,ij;
  if(A.sel==0 && B.sel==0) {
    if(A.n!=B.n || A.m!=B.m) on_error("matrix","diff. dimensions in A+B (1)");
    for(ij=0; ij<A.nm; ij++)  A.p[ij]+=B.p[ij];
  } else
  if(A.sel==0 && B.sel==1) {
    if(A.n!=B.n || A.m!=B.n) on_error("matrix","diff. dimension in A+B (2)");
    for(i=ij=0; i<A.n; i++)  for(j=0; j<A.m; j++,ij++)
      if(i==j)  A.p[ij]+=B.p[i];
  } else
  if(A.sel==1 && B.sel==1) {
    if(A.n!=B.n) on_error("matrix","different dimension in A+B (3)");
    for(i=0; i<A.n; i++)  A.p[i]+=B.p[i];
  } else
  if(A.sel==1 && B.sel==0)  return B+A;
  return A;
}
matrix operator-(matrix A)
{
  int ij;
  for(ij=0; ij<A.nm; ij++)  A.p[ij]=-A.p[ij];
  return A;
}
matrix operator-(complex x, matrix A)
{
  int i,j,ij;
  if(A.sel==0) {
    if(A.n!=A.m)  on_error("matrix","A is not square in x-A");
    for(i=ij=0; i<A.n; i++)  for(j=0; j<A.m; j++,ij++)
      if(i==j)  A.p[ij]=x-A.p[ij];  else  A.p[ij]=-A.p[ij];
  } else
  if(A.sel==1)  for(i=0; i<A.n; i++)  A.p[i]=x-A.p[i];
  return A;
}
matrix operator-(matrix A, const matrix &B)
{
  int i,j,ij;
  if(A.sel==0 && B.sel==0) {
    if(A.n!=B.n || A.m!=B.m) on_error("matrix","diff. dimensions in A-B (1)");
    for(ij=0; ij<A.nm; ij++)  A.p[ij]-=B.p[ij];
  } else
  if(A.sel==0 && B.sel==1) {
    if(A.n!=B.n || A.m!=B.n) on_error("matrix","diff. dimension in A-B (2)");
    for(i=ij=0; i<A.n; i++)  for(j=0; j<A.m; j++,ij++)
      if(i==j)  A.p[ij]-=B.p[i];
  } else
  if(A.sel==1 && B.sel==1) {
    if(A.n!=B.n) on_error("matrix","different dimension in A-B (3)");
    for(i=0; i<A.n; i++)  A.p[i]-=B.p[i];
  } else
  if(A.sel==1 && B.sel==0)  return -(B-A);
  return A;
}
matrix operator*(complex x, matrix A)
{
  int ij;
  for(ij=0; ij<A.nm; ij++)  A.p[ij]=x*A.p[ij];
  return A;
}
matrix operator*(numero x, matrix A)
{
  int ij;
  for(ij=0; ij<A.nm; ij++)  A.p[ij]=x*A.p[ij];
  return A;
}
matrix operator*(const matrix &A, const matrix &B)
{
  int i,j,k,ij;  complex val;  matrix C;
  complex *Cp,*Ap0,*Ap,*Bp,*Bp0;
  if(A.sel==0 && B.sel==0) {
    if(A.m!=B.n)  on_error("matrix","incompatible dimensions in A*B");
    C.alloc(A.n,B.m);
    for(i=0, Cp=C.p, Ap0=A.p; i<A.n; i++, Ap0+=A.m)
    for(j=0, Bp0=B.p; j<B.m; j++, Cp++, Bp0++) {
      for(k=0, val=0,  Ap=Ap0, Bp=Bp0;  k<A.m; k++, Ap++, Bp+=B.m)
        val+=(*Ap)*(*Bp);
      *Cp=val;
  } } else
  if(A.sel==0 && B.sel==1) {
    if(A.m!=B.n)  on_error("matrix","incompatible dimensions in A*B");
    C.alloc(A.n,A.m);
    for(i=ij=0; i<A.n; i++) for(j=0; j<A.m; j++,ij++) C.p[ij]=A.p[ij]*B.p[j];
  } else
  if(A.sel==1 && B.sel==0) {
    if(A.n!=B.n)  on_error("matrix","incompatible dimensions in A*B");
    C.alloc(B.n,B.m);
    for(i=ij=0; i<B.n; i++) for(j=0; j<B.m; j++,ij++) C.p[ij]=A.p[i]*B.p[ij];
  } else
  if(A.sel==1 && B.sel==1) {
    if(A.n!=B.n)  on_error("matrix","incompatible dimensions in A*B");
    C.alloc_diagonal(A.n);
    for(i=0; i<A.n; i++) C.p[i]=A.p[i]*B.p[i];
  }
  return C;
}
matrix transpose(matrix A)
{
  if(A.sel==1)  return A;
  matrix C(A.m,A.n);
  int i,j,ij;
  for(ij=0; ij<A.nm; ij++) {
    j=ij%A.m;  i=ij/A.m;
    C.p[j*C.m+i]=A.p[ij];
  }
  return C;
}
matrix conj(matrix A)
{
  matrix C(A.m,A.n);
  int i,j,ij;
  if(A.sel==0)
    for(ij=0; ij<A.nm; ij++) {
      j=ij%A.m;  i=ij/A.m;
      C.p[j*C.m+i]=conj(A.p[ij]);
    }
  else if(A.sel==1) for(ij=0; ij<A.nm; ij++) C.p[ij]=conj(A.p[ij]);
  return C;
}
matrix neg_off_diag(matrix A)
{
  matrix C(A.m,A.n);
  int i,j,ij,l=A.m/2;
  if(A.m!=A.n)  on_error("matrix::neg_off_diag", "n!=m");
  if(A.m%2)     on_error("matrix::neg_off_diag", "odd dimension");
  for(i=ij=0; i<A.n; i++) for(j=0; j<A.m; j++,ij++)
    if((i<l && j>=l) || (i>=l && j<l))  C.p[ij]=-A.p[ij];
    else                                C.p[ij]= A.p[ij];
  return C;
}
matrix &matrix::operator++(void)        
{
  int i,j,ij;
  if(sel==0)
    for(ij=0; ij<nm; ij++) {
      j=ij%m;  i=ij/m;
      p[j*m+i]=conj(p[ij]);
    }
  else if(sel==1) for(ij=0; ij<nm; ij++) p[ij]=conj(p[ij]);
  return *this;
}
matrix &matrix::operator--(void)        
{
  int ij;  if(sel==1) {for(ij=0; ij<nm; ij++) p[ij]=1/p[ij];  return *this;}
  matrix B;
  gaussj(B,n,0);
  return *this;
}
complex matrix::inv(void)               
{
  int ij;
  complex val;
  if(sel==1) {
    val=0;
    for(ij=0; ij<nm; ij++) {val+=log(p[ij]);  p[ij]=1/p[ij];}
    return val;
  }
  matrix B;
  return gaussj(B,n,0);
}
matrix &matrix::operator=(const matrix &A)
{
  int ij;
  if(A.nm)  alloc(A.n,A.m);  sel=A.sel;
  for(ij=0; ij<nm; ij++)  p[ij]=A.p[ij];
  return *this;
}
matrix &matrix::operator=(complex x)
{
  int i,j,ij;
  if(sel==0)
    for(i=ij=0; i<n; i++) for(j=0; j<m; j++,ij++)
      if(i==j) p[ij]=x; else  p[ij]=0;
  else if(sel==1) for(i=0; i<n; i++) p[i]=x;
  return *this;
}
matrix &matrix::operator=(numero x)
{
  int i,j,ij;
  if(sel==0)
    for(i=ij=0; i<n; i++) for(j=0; j<m; j++,ij++)
      if(i==j) p[ij]=x; else  p[ij]=0;
  else if(sel==1) for(i=0; i<n; i++) p[i]=x;
  return *this;
}
matrix operator%(matrix A, matrix B)
{
  int i,j,ij;
  if(A.sel==1 || B.sel==1) on_error("matrix","% not allowed with diag. mat.");
  if(A.nm==B.nm)                                               
    for(ij=0; ij<A.nm; ij++)  A.p[ij]=A.p[ij]*B.p[ij];  else
  if(A.nm==B.n) {                                              
    for(i=ij=0; i<B.n; i++)
    for(j=0; j<B.m; j++,ij++)  B.p[ij]=A.p[i]*B.p[ij];
    return B;
  }  else
  if(A.m==B.nm) {                                              
    for(i=ij=0; i<A.n; i++)
    for(j=0; j<A.m; j++,ij++)  A.p[ij]=A.p[ij]*B.p[j];
  }  else
  on_error("matrix","incompatible dimensions in A%B");
  return A;
}
matrix operator*(numero *r, matrix A)
{
  int ij;  complex *p;  p=A.p;
  for(ij=0; ij<A.nm; ij++, r++, p++)  *p=(*r)*(*p);
  return A;
}
matrix operator%(matrix A, numero *r)
{
  complex *p;  p=A.p;
  numero *rr;
  int i,j;
  for(i=0; i<A.n; i++) {
    rr=r;
    for(j=0; j<A.m; j++,p++,rr++)  *p=(*rr)*(*p);
  }
  return A;
}
matrix operator%(numero *r, matrix A)
{
  complex *p;  p=A.p;
  int i,j;
  for(i=0; i<A.n; i++, r++)
  for(j=0; j<A.m; j++,p++)  *p=(*r)*(*p);
  return A;
}
complex prod(matrix A, matrix B)
{
  complex val=0;
  int i;
  if(A.m!=1 || B.m!=1 || A.n!=B.n) on_error("matrix","A*B is not a scalar");
  for(i=0; i<A.n; i++)  val+=conj(A.p[i])*B.p[i];
  return val;
}
complex trace(const matrix &A)
{
  if(A.n!=A.m)  on_error("matrix","A is not square in trace(A)");
  int ij;  complex x(0,0);
  for(ij=0; ij<A.nm; ij+=A.n)  x+=A.p[ij];
  return x;
}
complex sum(const matrix &A)
{
  if(A.m!=1)  on_error("matrix","A is not (mx1) in sum(A)");
  int i;  complex x=0;
  for(i=0; i<A.n; i++)  x+=A.p[i];
  return x;
}
numero mod2(const matrix &A)
{
  int ij;  numero x=0;
  for(ij=0; ij<A.nm; ij++)  x+=mod2(A.p[ij]);
  return x;
}
complex det(const matrix &A)  {return exp(lndet(A));}
complex lndet(const matrix &A)
{
  matrix   B,C;   C=A;
  return C.gaussj(B,C.n,0,1);
}
void matrix::print(numero val)
{
  int ij;
  for(ij=0; ij<nm; ij++)  if(mod(p[ij])>=val)  printf("%3d %3d %15g %15g\n",
                                 ij/m, ij%m, real(p[ij]), imag(p[ij]));
}
complex matrix::gaussj(const matrix &Mb, int nn, int mm, int opt)
{
          if(n!=m || n<nn)  on_error("matrix","A not square in gaussj(A)");
  if(mm)  if(m!=Mb.n || Mb.m<mm)
            on_error("matrix","incompatible dimensions in gaussj");
  int  i,icol,irow,j,k,l,ll, *indxc,*indxr,*ipiv,*iA,*iB;
  numero  big,val;
  complex  dum,pivinv, *pirow,*picol,*pp,*ppp,*pa,*pb,*ppb, det=0;
  ipiv=new int [nn];  indxc=new int [nn];  indxr=new int [nn];
  iA=new int [nn];  iB=new int [nn];
  for(i=0; i<nn; i++)  {ipiv[i]=0;  iA[i]=i*m;  iB[i]=i*Mb.m;}
  for(i=0; i<nn; i++) {                    
    big=0;
    for(j=0,pp=p; j<nn; j++,pp+=m) {       
      if(ipiv[j]!=1)
      for(k=0,ppp=pp; k<nn; k++,ppp++) {
        if(!ipiv[k]) {
          val=mod(*ppp);
          if(val>=big) {
            big=val;
            irow=j;
            icol=k;
        } }  else  if(ipiv[k]>1)
                     on_error("matrix","1) singular matrix in gaussj");
    } }
    ipiv[icol]=ipiv[icol]+1;
    if(irow!=icol) {
      det=det+complex(0,pi);
      pirow=p+iA[irow];  picol=p+iA[icol];
      for(l=0; l<nn; l++,pirow++,picol++) {
        dum=*pirow;  *pirow=*picol;  *picol=dum;
      }
      if(mm)  {
        pirow=Mb.p+iB[irow];  picol=Mb.p+iB[icol];
        for(l=0; l<mm; l++,pirow++,picol++) {
          dum=*pirow;  *pirow=*picol;  *picol=dum;
    } } }
    indxr[i]=irow;
    indxc[i]=icol;   pp=p+iA[icol]+icol;
    if(real(*pp)==0 && imag(*pp)==0) {
      if(opt)  return -infinity;
      else     on_error("matrix","2) singular matrix in gaussj");
    }
    det=det+log(*pp);
    pivinv=1/(*pp);
    *pp=1;
    for(l=0, pp=p+iA[icol]; l<nn; l++,pp++)  *pp=(*pp)*pivinv;
    if(mm)  for(l=0, pp=Mb.p+iB[icol]; l<mm; l++,pp++)  *pp=(*pp)*pivinv;
    for(ll=0,pp=p,ppb=Mb.p; ll<nn; ll++,pp+=m,ppb+=mm)  if(ll!=icol)  {
        dum=*(pp+icol);
        *(pp+icol)=0;  pa=p+iA[icol];  pb=Mb.p+iB[icol];
        for(l=0,ppp=pp; l<nn; l++,ppp++,pa++)  *ppp=*ppp-(*pa)*dum;
        if(mm)  for(l=0,ppp=ppb; l<mm; l++,ppp++,pb++)  *ppp=*ppp-(*pb)*dum;
    }
  }
  for(l=nn-1; l>=0; l--)  if(indxr[l]!=indxc[l])
  for(k=0,pp=p+indxr[l],ppp=p+indxc[l]; k<nn; k++,pp+=m,ppp+=m) {
      dum=*pp;  *pp=*ppp;  *ppp=dum;
  }
  delete [] ipiv;  delete [] indxc;  delete [] indxr;
  delete [] iA;  delete [] iB;
  val=imag(det);  val-=int(val/(2*pi))*2*pi;            
  while(val<-2*pi) val+=2*pi;  while(val>2*pi) val-=2*pi;
  return complex(real(det), val);
}
void escribe_matrix(matrix &a)
{
  int i,j,ij;
  for(i=ij=0; i<a.n; i++)
  for(j=0; j<a.m; j++, ij++) {
    printf("  (%2d:%2d) %14.6g", i, j, mod(a.a(i,j)));
    if((ij+1)%3==0)  printf("\n");
  }
  printf("\n");
}
#endif  
#ifdef jga_recurs         
#else                     
#define jga_recurs 1      
                          
                          
class recurs {
private:
  complex bb(int p);                     
  complex Ij(int j, int p, int q);       
  complex *I, *a, *b;
  int     *iI, *ib, n, n2, nn,
          init_flag;                     
public:
  int init(int n_, complex *mu);
  void free(void);
  complex eval(complex E, int n_, int termination);
  complex eval(numero  E, int n_, int termination);
  complex eval(complex E, int n_, complex *mu, int termination);
  complex eval(numero  E, int n_, complex *mu, int termination);
  recurs(void) {init_flag=0;}
};
complex recurs::bb(int p)
{
  if(p<1)  return 0;
  if(!ib[p]) {
    b[p]=sqrt(Ij(2,p-1,p-1)-sqr(Ij(1,p-1,p-1))-sqr(bb(p-1)));
    ib[p]=1;
  }
  return b[p];
}
complex recurs::Ij(int j, int p, int q)
{
  if(p<0 || q<0)  return 0;
  int k=j*n2+p*n+q;  complex bb_;
  if(!iI[k]) {
    if(p>0)  bb_=bb(p);  else  bb_=bb(q);
    if(ABS(real(bb_))<infinitesimal &&
       ABS(imag(bb_))<infinitesimal)   I[k]=0; else
    if(p>0)  I[k]=( Ij(j+1,p-1,q)
                   -Ij(1,p-1,p-1)*Ij(j,p-1,q)
                   -bb(p-1)*Ij(j,p-2,q) )/bb_;
    else     I[k]=( Ij(j+1,p,q-1)
                   -Ij(1,q-1,q-1)*Ij(j,p,q-1)
                   -bb(q-1)*Ij(j,p,q-2) )/bb_;
    iI[k]=1;
  }
  return I[k];
}
int recurs::init(int n_, complex *mu)
{
  if(n_>=1)
  if(mod(mu[0])<1e-10*mod(mu[1]) || mod(mu[0])<1e-10) {   
    mu++;
    init(n_-1,mu);
    return 0;
  }
  free();  init_flag=1;  nn=n_;
  if(nn<1) {
    n=0;
    b=new complex [1+1000];   b[0]=mu[0];
    return 0;
  }
  n =(nn-1)/2 + 1;
  n2=sqr(n);
  int i,n3=(2*n+1)*n2+n2+n;
  a=new complex [n+2];
  b=new complex [n+2];
  ib=new int [n+2];
  I=new complex [n3];
  iI=new int [n3];
  for(i=0; i<n; i++)  ib[i]=0;
  for(i=0; i<n3; i++)  iI[i]=0;
  for(i=0; i<=nn; i++) {
    I[i*n2]=mu[i]/mu[0];
    iI[i*n2]=1;
  }
  for(i=0; i<n; i++)  a[i]=Ij(1,i,i);
  for(i=1; i<n; i++)  b[i]=bb(i);
  for(i=1; i<n; i++)  b[i]=sqr(b[i]);
                      b[0]=mu[0];
  delete [] I;   delete [] iI;   delete [] ib;
  return 0;
}
void recurs::free(void)
{
  if(init_flag)  {if(n>0) delete [] a;  delete [] b;}
}
complex recurs::eval(complex E, int n_, int termination)
{
  if(n_>nn)  n_=nn;
  int  nt=(n_-1)/2 + 1;
  if(n_==0)  return b[0]/E;
  complex t;                      
  int i;                          
  if(!termination || nt==1)  t=0;
  else  t=(E-a[nt-1]-sqrt(sqr(E-a[nt-1]) - 4*b[nt-1]))/2.0;
  for(i=nt-1; i>=0; i--)  t=b[i]/(E-a[i]-t);
  return  t;
}
complex recurs::eval(complex E, int n_, complex *mu, int termination)
{
  complex val;
  init(n_,mu);
  val=eval(E,n_,termination);
  free();
  return val;
}
complex recurs::eval(numero E, int n_, int termination)
  {return  eval(complex(E,0),n_,termination);}
complex recurs::eval(numero E, int n_, complex *mu, int termination)
  {return  eval(complex(E,0),n_,mu,termination);}
#endif  
#ifdef jga_light          
#else                     
#define jga_light 1       
                          
                          
void light_polarization_mu(numero th, numero fi, numero al, numero del,
                           complex *eps)
{
  numero cal=cos(al), cth=cos(th), sal=sin(al), sth=sin(th);
  eps[0]=( cal*cth + i_c*sal*exp(i_c*del))*exp( i_c*fi)/sqrt(2);  
  eps[1]= -cal*sth;                                               
  eps[2]=(-cal*cth + i_c*sal*exp(i_c*del))*exp(-i_c*fi)/sqrt(2);  
}
void light_polarization_xyz(complex *eps,
                            complex &ex, complex &ey, complex &ez)
{
  ex=    -(eps[2]-eps[0])/sqrt(2);
  ey=-i_c*(eps[2]+eps[0])/sqrt(2);
  ez=      eps[1];
}
void light_polarization_xyz(numero th, numero fi, numero al, numero del,
                            complex *EE)
{
  complex eps[3], ex,ey,ez;
  light_polarization_mu(th,fi,al,del, eps);
  light_polarization_xyz(eps, ex, ey, ez);
  EE[0]=ex;  EE[1]=ey;  EE[2]=ez;
}
void light_E_field(numero th, numero fi, numero alpha, numero delta,
                   complex &ex, complex &ey, complex &ez)
{
  ex= cos(alpha)*cos(th)*cos(fi)-sin(alpha)*sin(fi)*exp(i_c*delta);
  ey= cos(alpha)*cos(th)*sin(fi)+sin(alpha)*cos(fi)*exp(i_c*delta);
  ez=-cos(alpha)*sin(th);
}
void light_A_potential(complex w,
                       numero th, numero fi, numero alpha, numero delta,
                       complex &ax, complex &ay, complex &az)
{
  complex ik=i_c*w/c_au;
  light_E_field(th,fi,alpha,delta,ax,ay,az);
  ax=ax/ik;  ay=ay/ik;  az=az/ik;
}
#endif  
#ifndef jga_Salvat        
#define jga_Salvat 1      
                          
                          
class Salvat {
public:
numero *r, *rr, *P, *Q;
int    ngp, nmax, ndim;
spline PP,QQ;
numero phase_last, E_last, kk_last;                  
int    k_last;                                       
void vint(numero *R_, numero *RV_, int np);
void vint(numero *r_in, numero *V_in, int n_in, int atZ);
void dbound(numero &E, numero err, int N, int K);
void dfree(numero E, numero &PHASE, int k);
void printwf(FILE *fout);
void free(void);
void wf(numero r_, numero &f, numero &g);
void phases(numero *r_in, numero *V_in, int n_in,    
            int atZ,                                 
            numero Ei, numero Ef, int nE,            
            int rel, int lmax, numero *ph);
void phases(numero *r_in, numero *V_in, int n_in,    
            int atZ, numero E,                       
            int rel, int lmax, numero *ph)
  {phases(r_in, V_in, n_in, atZ, E,E,1, rel, lmax, ph);}
numero eps;
Salvat(void) {eps=1.0e-15;  ndim=0;}
numero *RV, *VA, *VB, *VC, *VD;
numero *Y, *A, *B, *C, *D;
int    NZMAX, NSTEP,NCHS, NSUM;
numero RV0, RV1, RV2, RV3;
numero PPI,QQI,PF,QF,RA,RB,RLN;
numero P0,Q0,P1,Q1,CA[61],CB[61],R0,R1;
numero vacuum(int l, numero kr, numero ps)
  {return besselj(l,kr)*cos(ps) - bessely(l,kr)*sin(ps);}
int    DOUTW(numero E, int K, int NR, int &NCERO, int IOTP);
void   DINW(numero E, int K, int IOTP);
void   DIR(numero E, numero AK);
int    DIR0(numero E, numero AK);
void   SPLINE(numero *X, numero *Y, numero S1, numero SN, int &N);
numero INTEG(numero *X, numero XL, numero XU, int N);
int    FINDI(numero *X, numero XC, int N);
numero DMAX1(numero a, numero b) {if(a>b) return a;  else  return b;}
numero DMIN1(numero a, numero b) {if(a<b) return a;  else  return b;}
numero DMAX1(numero a, numero b, numero c, numero d)
       {return DMAX1(DMAX1(a,b),DMAX1(c,d));}
};
void Salvat::free(void) {
  if(ndim) {
    delete [] r;  delete [] rr;  delete [] P;  delete [] Q;
    delete [] RV;  delete [] VA;  delete [] VB;  delete [] VC;  delete [] VD;
    delete [] Y;  delete [] A;  delete [] B;  delete [] C;  delete [] D;
  }
  PP.free();  QQ.free();  ndim=0;
}
void Salvat::vint(numero *R_, numero *RV_, int np)
{
  free();  ndim=np+1;
  r=new numero [ndim];  rr=new numero [ndim];
  P=new numero [ndim];  Q=new numero [ndim];
  RV=new numero [ndim];  VA=new numero [ndim];  VB=new numero [ndim];
  VC=new numero [ndim];  VD=new numero [ndim];
  Y=new numero [ndim];  A=new numero [ndim];  B=new numero [ndim];
  C=new numero [ndim];  D=new numero [ndim];
  int IO=0,I=0,K=0,J;                          
gt1:    I++;  K++;                             
        P[K]=R_[I-1];  Q[K]=RV_[I-1];          
        if(I==np) goto gt2;
        if(R_[I-1]<R_[I]-1.0e-12) goto gt1;
gt2:    SPLINE(P,Q,0,0,K);
        K--;
        for(J=1; J<=K; J++) {
          IO++;
          rr[IO-1]=r[IO]=P[J];
          RV[IO]=Q[J];
          VA[IO]=A[J];  VB[IO]=B[J];
          VC[IO]=C[J];  VD[IO]=D[J];
        }
  if(I<np) {K=0; goto gt1;}
  ngp=IO+1;
  rr[ngp-1]=r[ngp]=P[K+1];
  RV[ngp]=Q[K+1];
}
void Salvat::wf(numero r_, numero &f, numero &g)
{
  if(ndim==0)  on_error("Salvat::wf(r,f,g)","no wf is not available");
  numero kr;
  int l,lp;
  if(r_<=PP.b) {f=PP.val(r_);  g=QQ.val(r_);}  else
  if(E_last>0) {
    kr=kk_last*r_;
    if(k_last<0) l=-k_last-1; else l=k_last;
    if(k_last<0) lp=l+1; else lp=l-1;
    f = r_ * vacuum(l,kr,phase_last);
    g = -r_ * vacuum(lp,kr,phase_last) * E_last/(kk_last*c_au);
  }  else  on_error("Salvat::wf(r,f,g)","r out of range for bound state");
}
void Salvat::dbound(numero &E, numero err, int N, int K)
{
      numero FL1,ESUP,EMIN,EMAX,EKIN,PO,QO,FACT,RLAST,SUM,DE,EP,EO;
      int NR,ICMIN,ICMAX,L,I,J,IOTP,NCERO;
      k_last=K;   phase_last=0;
      if(N<1||K==0||K<-N||K>N-1)
        on_error("Salvat::dbound","n or k out of range");
      if(K<0) L=-K-1; else L=K;
      NR=N-L-1;
      err=ABS(err);
      if(E>-0.1) E=-0.1;
      FL1=0.5*L*(L+1);
      EMIN=1;
      for(I=2; I<=ngp; I++) EMIN=DMIN1(EMIN,(RV[I]+FL1/r[I])/r[I]);
      if(EMIN>=0)
        on_error("Salvat::dbound", "error 1: EMIN>0: use denser grid");
      ESUP=(RV[ngp]+FL1/r[ngp])/r[ngp];
      if(ESUP>0) ESUP=0;
      if(EMIN<-c_au*c_au) EMIN=-c_au*c_au;
      if(E>ESUP || E<EMIN) E=0.5*(ESUP+EMIN);
      EMAX=ESUP;
      ICMIN=0;
      ICMAX=0;
gt2:  if(E>-1.0e-16)
        on_error("Salvat::dbound", "error 2: bound state does not exist");
      for(J=2; J<=ngp; J++) {
        IOTP=ngp+2-J;
        EKIN=E-(RV[IOTP]+FL1/r[IOTP])/r[IOTP];
        if(EKIN>0) goto gt4;
      }
gt4:  IOTP++;
      if(IOTP>ngp-1)  on_error("Salvat::dbound", "error 3: use denser grid");
      DOUTW(E,K,NR,NCERO,IOTP);
      if(NZMAX>1 && NCERO<=NR)
        on_error("Salvat::dbound", "error 4: use denser grid");
      if(NCERO>NR) {
        if(ICMIN==0) EMIN=EMIN-(EMAX-E);
        EMAX=E;
        ICMAX=1;
        E=0.5*(EMIN+E);
        if(EMAX-EMIN<err*ABS(EMIN))
          on_error("Salvat::dbound", "error 5: E out of range");
        goto gt2;
      }
      if(NCERO<NR) {
        ICMIN=1;
        EMIN=E;
        E=0.5*(E+EMAX);
        if(EMAX-EMIN<err*ABS(EMIN))
          on_error("Salvat::dbound", "error 5: E out of range");
        goto gt2;
      }
      PO=P[IOTP];
      QO=Q[IOTP];
      DINW(E,K,IOTP);
      FACT=PO/P[IOTP];
      for(I=IOTP; I<=ngp; I++) {
        P[I]=P[I]*FACT;
        Q[I]=Q[I]*FACT;
      }
      QQI=Q[IOTP];
      for(I=1; I<=nmax; I++) Y[I]=sqr(P[I])+sqr(Q[I]);
      RLAST=r[nmax];
      SPLINE(r,Y,0,0,nmax);
      SUM=INTEG(r,0,RLAST,nmax);
      if(SUM<1.0e-15) SUM=1;
      DE=c_au*PO*(QQI-QO)/SUM;
      EP=E+DE;
      if(DE<0) {ICMAX=1; EMAX=E;}
      if(DE>0) {ICMIN=1; EMIN=E;}
      if(ICMIN==0 && EP<EMIN) {
        EMIN=1.1*EMIN;
        if(EP<EMIN) EP=0.5*(E+EMIN);
      }
      if(ICMIN==1 && EP<EMIN) EP=0.5*(E+EMIN);
      if(ICMAX==1 && EP>EMAX) EP=0.5*(E+EMAX);
      if(EP>ESUP) EP=0.5*(ESUP+E);
      if(EP>ESUP && ABS(E-ESUP)<err*ABS(ESUP))
        on_error("Salvat::dbound", "error 5: E out of range");
      EO=E;
      E=EP;
      if(DMIN1(ABS(DE),ABS(E-EO))>ABS(E*err)) goto gt2;
      FACT=1/sqrt(SUM);
      PP.alloc(nmax);  QQ.alloc(nmax);
      for(I=0; I<nmax; I++) {
        P[I]=P[I+1]*FACT; Q[I]=-Q[I+1]*FACT;
        PP.put(I,rr[I],P[I]);  QQ.put(I,rr[I],Q[I]);
      }
      PP.init();  QQ.init();
      for(I=nmax; I<ndim; I++)  P[I]=Q[I]=0;
      if(ABS(P[ngp])>1.0e-3*ABS(P[IOTP]))
        on_error("Salvat::dbound", "error 3: use larger radii");
      E_last=E;
}
void Salvat::dfree(numero E, numero &PHASE, int K)
{
      int L,I,IL,NCERO;
      numero FL1,T,RF,PF_,QF_,RK,BJL,BJLP,BNL,BNLP,TT,VF,PPF,RATIO,FACT;
      if(K==0)  on_error("Salvat::dfree", "K=0");  k_last=K;
      if(E<=0)  on_error("Salvat::dfree", "E<0");  E_last=E;
      if(K<0) L=-K-1; else L=K;
      FL1=0.5*L*(L+1);
      T=1.0e-6*r[ngp]*ABS(E-FL1/(r[ngp]*r[ngp]));
      if(ABS(RV[ngp])>T)
        on_error("Salvat::dfree", "error 7: extend grid to larger radii");
      nmax=ngp;
      for(I=4; I<=ngp; I++) {
        IL=nmax-1;
        T=1-10*r[IL]*ABS(E-FL1/(r[IL]*r[IL]));
        if(ABS(RV[IL])>T) goto gt2;
        nmax=IL;
      }
gt2:
      DOUTW(E,K,1,NCERO,nmax);
      RF=r[nmax];
      PF_=P[nmax];
      QF_=Q[nmax];
      RK=sqrt(E*(E+2*sqr(c_au)))/c_au;
      T=RK*RF;
      BJL=besselj(L,T);           BNL=bessely(L,T);
      besseljyp(L,T,BJLP,BNLP);
      IL=nmax-1;
      VF=VA[IL]/RF+VB[IL]+RF*(VC[IL]+RF*VD[IL]);
      PPF=-K*PF_/RF-(E-VF+2*sqr(c_au))*QF_/c_au;
      RATIO=PPF/PF_-1/RF;
      PHASE=atan2(RK*BJLP-RATIO*BJL,RK*BNLP-RATIO*BNL);
      TT=ABS(PHASE);
      if(TT>0.5*pi) PHASE=PHASE*(1-pi/TT);
      kk_last=sqrt(E*(E+2*c_au*c_au))/c_au;
      FACT=r[nmax]*vacuum(L,kk_last*r[nmax],PHASE)/P[nmax];
      PP.alloc(nmax);  QQ.alloc(nmax);
      for(I=0; I<nmax; I++) {                  
        P[I]=P[I+1]*FACT; Q[I]=-Q[I+1]*FACT;
        PP.put(I,rr[I],P[I]);  QQ.put(I,rr[I],Q[I]);
      }
      PP.init();  QQ.init();
      for(I=nmax; I<ndim; I++)  P[I]=Q[I]=0;
      phase_last=PHASE;
}
int Salvat::DOUTW(numero E, int K, int NR, int &NCERO, int IOTP)
{
  numero AK,FACT;
  int N1,I,J;
      NCERO=0;
      NZMAX=0;
      AK=K;
      N1=IOTP-1;
      P[1]=0;
      Q[1]=0;
      for(I=1; I<=N1; I++) {
        RA=r[I];
        RB=r[I+1];
        RV0=VA[I];
        RV1=VB[I];
        RV2=VC[I];
        RV3=VD[I];
        PPI=P[I];
        QQI=Q[I];
        DIR(E,AK);
        NCERO=NCERO+NCHS;
        if(NCHS>NZMAX) NZMAX=NCHS;
        if(NCERO>NR && E<0) return 1;
        P[I+1]=PF;
        Q[I+1]=QF;
        if(I>1 && RLN>0) {
          FACT=exp(-RLN);
          for(J=1; J<=I; J++) {
            P[J]=P[J]*FACT;
            Q[J]=Q[J]*FACT;
      } } }
      return 1;
}
void Salvat::DINW(numero E, int K, int IOTP)
{
      numero AK,AL,RN,RVN,RVNP,CMU,CRAT,FACT, TRINF=5625.0;
      int L,N,N1,I,J,I1,M;
      if(K<0) L=-K-1; else L=K;
      AK=K;
      AL=L;
      N=ngp;
      FACT=(E+2*c_au*c_au)/(c_au*c_au);
gt1:  N1=N-1;
      RN=r[N];
      RVN=VA[N1]+RN*(VB[N1]+RN*(VC[N1]+RN*VD[N1]));
      RVNP=VB[N1]+RN*(2*VC[N1]+RN*3*VD[N1]);
      CMU=FACT*RN*(RVN-E*RN)+AL*(AL+1);
      if(CMU<=0)  on_error("Salvat::DINW", "error 6: RV[ngp]<<0 or E>0");
      if(CMU<TRINF || N==IOTP+1) {
        CRAT=(0.5-sqrt(CMU))/RN-0.25*FACT*(RVN+RN*RVNP-2*RN*E)/CMU;
        P[N]=1;
        Q[N]=-c_au*(CRAT+AK/RN)/(E+2*c_au*c_au);
        nmax=N;
      } else {
        P[N]=0;
        Q[N]=0;
        N--;
        goto gt1;
      }
      N1=N-IOTP;
      for(J=1; J<=N1; J++) {
        I=N-J;
        I1=I+1;
        RA=r[I1];
        RB=r[I];
        RV0=VA[I];
        RV1=VB[I];
        RV2=VC[I];
        RV3=VD[I];
        PPI=P[I1];
        QQI=Q[I1];
        DIR(E,AK);
        P[I]=PF;
        Q[I]=QF;
        if(RLN>0) {
          FACT=exp(-RLN);
	  for(M=I1; M<=N; M++) {
            P[M]=P[M]*FACT;
            Q[M]=Q[M]*FACT;
}     } } }
void Salvat::DIR(numero E, numero AK)
{
      numero H, DIRECT,TST;
      int K,IOUT;
      NCHS=0;
      RLN=0;
      H=RB-RA;
      if(H<0) DIRECT=-1; else DIRECT=1;
      K=-2;
      NSTEP=0;
      R1=RA;
      P1=PPI;
      Q1=QQI;
gtD1: R0=R1;
      P0=P1;
      Q0=Q1;
gtD2: IOUT=0;
      R1=R0+H;
      if(DIRECT*(RB-R1)<DIRECT*0.1*H) {
        R1=RB;
        H=RB-R0;
        IOUT=1;
      }
      DIR0(E,AK);
      K++;
      if(NSUM>15) goto gtD3;
      if(K<0) goto gtD4;
      H+=H;
      K=0;
      goto gtD4;
gtD3: if(NSUM<60) goto gtD4;
      H=0.5*H;
      K=-4;
      goto gtD2;
gtD4: NSTEP++;
      TST=ABS(P1);
      if(TST>100) {
        RLN=RLN+log(TST);
        P1=P1/TST;
        Q1=Q1/TST;
      }
      if(P0*P1<0 && R0>0) NCHS++;
      if(IOUT==0) goto gtD1;
      PF=P1;
      QF=Q1;
}
int Salvat::DIR0(numero E, numero AK)
{
      int ISIG,I,K;
      numero H,H2,RVE,U0,U1,U2,U3,UT,UQ,UH,S,CAI,CBI,DS,QP1,PP1,TST,T1A,
             T1B,T1,T2,DS1,RHO,UB;
      numero OVER=1.0e15;
      ISIG=1;
      if(AK>0) ISIG=-1;
      H=R1-R0;
      H2=H*H;
      RVE=RV1-E;
      if(R0>1.0e-10) goto gt07;
      U0=RV0/c_au;
      U1=RVE*R1/c_au;
      U2=RV2*R1*R1/c_au;
      U3=RV3*R1*R1*R1/c_au;
      UT=U0+U1+U2+U3;
      UQ=UT-2*c_au*R1;
      UH=U1-2*c_au*R1;
      if(ABS(U0)<1.0e-10) goto gt02;
      S=sqrt(AK*AK-U0*U0);
      DS=S+S;
      CA[1]=1;
      CB[1]=(S+AK)/U0;
      CAI=U1*CA[1];
      CBI=UH*CB[1];
      CA[2]=(U0*CAI+(S+1-AK)*CBI)/(DS+1);
      CB[2]=(-(S+1+AK)*CAI+U0*CBI)/(DS+1);
      CAI=U1*CA[2]+U2*CA[1];
      CBI=UH*CB[2]+U2*CB[1];
      CA[3]=(U0*CAI+(S+2-AK)*CBI)/(2*(DS+2));
      CB[3]=(-(S+2+AK)*CAI+U0*CBI)/(2*(DS+2));
      P1=CA[1]+CA[2]+CA[3];
      PP1=S*CA[1]+(S+1)*CA[2]+(S+2)*CA[3];
      Q1=CB[1]+CB[2]+CB[3];
      QP1=S*CB[1]+(S+1)*CB[2]+(S+2)*CB[3];
      for(I=4; I<=60; I++) {
        K=I-1;
        CAI=U1*CA[K]+U2*CA[I-2]+U3*CA[I-3];
        CBI=UH*CB[K]+U2*CB[I-2]+U3*CB[I-3];
        CA[I]=(U0*CAI+(S+K-AK)*CBI)/(K*(DS+K));
        CB[I]=(-(S+K+AK)*CAI+U0*CBI)/(K*(DS+K));
        P1=P1+CA[I];
        PP1=PP1+(S+K)*CA[I];
        Q1=Q1+CB[I];
        QP1=QP1+(S+K)*CB[I];
        TST=DMAX1(ABS(P1),ABS(Q1),ABS(PP1),ABS(QP1));
        if(TST>OVER) {NSUM=100;  return 1;}
        T1A=ABS(R1*PP1+H*(AK*P1-UQ*Q1));
        T1B=ABS(R1*QP1-H*(AK*Q1-UT*P1));
        T1=DMAX1(T1A,T1B);
        T2=DMAX1(ABS(CA[I]),ABS(CB[I]));
        TST=eps*DMAX1(ABS(P1),ABS(Q1));
        if(T1<TST && T2<TST) goto gt06;
      }
      goto gt06;
gt02: if(ISIG<0) goto gt04;
      S=ABS(AK);
      DS1=S+S+1;
      CA[1]=1;
      CB[1]=-U1*CA[1]/DS1;
      CA[2]=0;
      CB[2]=-U2*CA[1]/(DS1+1);
      CA[3]=UH*CB[1]/2;
      CB[3]=-(U1*CA[3]+U3*CA[1])/(DS1+2);
      CA[4]=(UH*CB[2]+U2*CB[1])/3;
      CB[4]=-(U1*CA[4]+U2*CA[3])/(DS1+3);
      P1=CA[1]+CA[2]+CA[3]+CA[4];
      PP1=S*CA[1]+(S+1)*CA[2]+(S+2)*CA[3]+(S+3)*CA[4];
      Q1=CB[1]+CB[2]+CB[3]+CB[4];
      QP1=(S+1)*CB[1]+(S+2)*CB[2]+(S+3)*CB[3];
      for(I=5; I<=60; I++) {
        K=I-1;
        CA[I]=(UH*CB[I-2]+U2*CB[I-3]+U3*CB[I-4])/K;
        CB[I]=-(U1*CA[I]+U2*CA[K]+U3*CA[I-2])/(DS1+K);
        P1=P1+CA[I];
        PP1=PP1+(S+K)*CA[I];
        Q1=Q1+CB[I];
        QP1=QP1+(S+I)*CB[I];
        TST=DMAX1(ABS(P1),ABS(Q1),ABS(PP1),ABS(QP1));
        if(TST>OVER) {NSUM=100;  return 1;}
        T1A=ABS(R1*PP1+H*(AK*P1-UQ*Q1));
        T1B=ABS(R1*QP1-H*(AK*Q1-UT*P1));
        T1=DMAX1(T1A,T1B);
        T2=DMAX1(ABS(CA[I]),ABS(CB[I]));
        TST=eps*DMAX1(ABS(P1),ABS(Q1));
        if(T1<TST && T2<TST) goto gt06;
      }
      goto gt06;
gt04: S=ABS(AK)+1;
      DS1=S+ABS(AK);
      CB[1]=1;
      CA[1]=UH*CB[1]/DS1;
      CB[2]=0;
      CA[2]=U2*CB[1]/(DS1+1);
      CB[3]=-U1*CA[1]/2;
      CA[3]=(UH*CB[3]+U3*CB[1])/(DS1+2);
      CB[4]=-(U1*CA[2]+U2*CA[1])/3;
      CA[4]=(UH*CB[4]+U2*CB[3])/(DS1+3);
      P1=CA[1]+CA[2]+CA[3]+CA[4];
      PP1=S*CA[1]+(S+1)*CA[2]+(S+2)*CA[3]+(S+3)*CA[4];
      Q1=CB[1]+CB[2]+CB[3]+CB[4];
      QP1=(S-1)*CB[1]+S*CB[2]+(S+1)*CB[3];
      for(I=5; I<=60; I++) {
        K=I-1;
        CB[I]=-(U1*CA[I-2]+U2*CA[I-3]+U3*CA[I-4])/K;
        CA[I]=(UH*CB[I]+U2*CB[K]+U3*CB[I-2])/(DS1+K);
        P1=P1+CA[I];
        PP1=PP1+(S+K)*CA[I];
        Q1=Q1+CB[I];
        QP1=QP1+(S+K-1)*CB[I];
        TST=DMAX1(ABS(P1),ABS(Q1),ABS(PP1),ABS(QP1));
        if(TST>OVER) {NSUM=100;  return 1;}
        T1A=ABS(R1*PP1+H*(AK*P1-UQ*Q1));
        T1B=ABS(R1*QP1-H*(AK*Q1-UT*P1));
        T1=DMAX1(T1A,T1B);
        T2=DMAX1(ABS(CA[I]),ABS(CB[I]));
        TST=eps*DMAX1(ABS(P1),ABS(Q1));
        if(T1<TST && T2<TST) goto gt06;
      }
gt06: NSUM=K+1;
      Q1=Q1/P1;
      P1=1;
      return 1;
gt07:
      RHO=H/R0;
      U0=(RV0+R0*(RVE+R0*(RV2+R0*RV3)))/c_au;
      U1=(RVE+R0*(2*RV2+R0*3*RV3))*H/c_au;
      U2=(RV2+R0*3*RV3)*H2/c_au;
      U3=RV3*H*H2/c_au;
      UB=U0-2*c_au*R0;
      UH=U1-2*c_au*H;
      UT=U0+U1+U2+U3;
      UQ=UT-2*c_au*R1;
      CA[1]=P0;
      CB[1]=Q0;
      CA[2]=RHO*(-AK*CA[1]+UB*CB[1]);
      CB[2]=-RHO*(U0*CA[1]-AK*CB[1]);
      CA[3]=RHO*(-(1+AK)*CA[2]+UB*CB[2]+UH*CB[1])/2;
      CB[3]=-RHO*(U0*CA[2]+(1-AK)*CB[2]+U1*CA[1])/2;
      CAI=-(2+AK)*CA[3]+UB*CB[3]+UH*CB[2]+U2*CB[1];
      CBI=U0*CA[3]+(2-AK)*CB[3]+U1*CA[2]+U2*CA[1];
      CA[4]=RHO*CAI/3;
      CB[4]=-RHO*CBI/3;
      P1=CA[1]+CA[2]+CA[3]+CA[4];
      PP1=CA[2]+2*CA[3]+3*CA[4];
      Q1=CB[1]+CB[2]+CB[3]+CB[4];
      QP1=CB[2]+2*CB[3]+3*CB[4];
      for(I=5; I<=60; I++) {
        K=I-1;
        CAI=-(K-1+AK)*CA[K]+UB*CB[K]+UH*CB[I-2]+U2*CB[I-3]+U3*CB[I-4];
        CBI= U0*CA[K]+(K-1-AK)*CB[K]+U1*CA[I-2]+U2*CA[I-3]+U3*CA[I-4];
        CA[I]=RHO*CAI/K;
        CB[I]=-RHO*CBI/K;
        P1=P1+CA[I];
        PP1=PP1+K*CA[I];
        Q1=Q1+CB[I];
        QP1=QP1+K*CB[I];
        TST=DMAX1(ABS(P1),ABS(Q1),ABS(PP1),ABS(QP1));
        if(TST>OVER) {NSUM=100;  return 1;}
        T1A=ABS(R1*PP1+H*(AK*P1-UQ*Q1));
        T1B=ABS(R1*QP1-H*(AK*Q1-UT*P1));
        T1=DMAX1(T1A,T1B);
        T2=DMAX1(ABS(CA[I]),ABS(CB[I]));
        TST=eps*DMAX1(ABS(P1),ABS(Q1));
        if(T1<TST && T2<TST) goto gt010;
      }
gt010: NSUM=K+1;
      return 1;
}
void Salvat::SPLINE(numero *X, numero *Y_, numero S1, numero SN, int &N)
{
      int N1=N-1;
      int N2=N-2;
      int I,K;
      numero H,HI,R,SI,SI1;
      if(N<4)  on_error("Salvat::SPLINE", "less than 4 points");
      for(I=1; I<=N1; I++) {
        if(X[I+1]-X[I]<1.0e-10)
          on_error("Salvat::SPLINE", "X is not in increasing order");
        A[I]=X[I+1]-X[I];
        D[I]=(Y_[I+1]-Y_[I])/A[I];
      }
      for(I=1; I<=N2; I++) {
        B[I]=2.0*(A[I]+A[I+1]);
        K=N1-I+1;
        D[K]=6.0*(D[K]-D[K-1]);
      }
      D[2]=D[2]-A[1]*S1;
      D[N1]=D[N1]-A[N1]*SN;
      for(I=2; I<=N2; I++) {
        R=A[I]/B[I-1];
        B[I]=B[I]-R*A[I];
        D[I+1]=D[I+1]-R*D[I];
      }
      D[N1]=D[N1]/B[N2];
      for(I=2; I<=N2; I++) {
        K=N1-I+1;
        D[K]=(D[K]-A[K]*D[K+1])/B[K-1];
      }
      D[N]=SN;
      SI1=S1;
      for(I=1; I<=N1; I++) {
        SI=SI1;
        SI1=D[I+1];
        H=A[I];
        HI=1/H;
        A[I]= HI/6*(SI*sqr(X[I+1])*X[I+1]-SI1*sqr(X[I])*X[I])
              +HI*(Y_[I]*X[I+1]-Y_[I+1]*X[I])
              +H/6*(SI1*X[I]-SI*X[I+1]);
        B[I]= HI/2*(SI1*sqr(X[I])-SI*sqr(X[I+1]))
              +HI*(Y_[I+1]-Y_[I])+(H/6)*(SI-SI1);
        C[I]= HI/2*(SI*X[I+1]-SI1*X[I]);
        D[I]= HI/6*(SI1-SI);
}     }
numero Salvat::INTEG(numero *X, numero XL, numero XU, int N)
{
  numero SUM, SIGN, X1,X2, SUMP;
  int IL,IU,I;
      SIGN=1;
      if(XU<XL) {
        SUM=XL;
        XL=XU;
        XU=SUM;
        SIGN=-1;
      }
      if(XL<X[1] || XU>X[N])
        on_error("Salvat::INTEG", "integral limits out of range");
      SUM=0;
      IL=FINDI(X,XL,N);
      IU=FINDI(X,XU,N);
      if(IL==IU) {
        X1=XL;
        X2=XU;
        SUM= X2*(A[IL]+X2*((B[IL]/2)+X2*((C[IL]/3)+X2*D[IL]/4)))
            -X1*(A[IL]+X1*((B[IL]/2)+X1*((C[IL]/3)+X1*D[IL]/4)));
      } else {
        X1=XL;
        X2=X[IL+1];
        SUM= X2*(A[IL]+X2*((B[IL]/2)+X2*((C[IL]/3)+X2*D[IL]/4)))
            -X1*(A[IL]+X1*((B[IL]/2)+X1*((C[IL]/3)+X1*D[IL]/4)));
        IL=IL+1;
        for(I=IL; I<=IU; I++) {
          X1=X[I];
          X2=X[I+1];
          if(I==IU) X2=XU;
          SUMP= X2*(A[I]+X2*((B[I]/2)+X2*((C[I]/3)+X2*D[I]/4)))
               -X1*(A[I]+X1*((B[I]/2)+X1*((C[I]/3)+X1*D[I]/4)));
          SUM=SUM+SUMP;
      } }
      return SIGN*SUM;
}
int Salvat::FINDI(numero *X, numero XC, int N)
{
  int I=1, I1=N, IT;
  if(XC>X[N]) return I=N-1;
  if(XC<X[1]) return I=1;
  do {
    IT=(I+I1)/2;
    if(XC>X[IT])  I=IT;
    if(XC<=X[IT]) I1=IT;
  } while(I1-I>1);
  return I;
}
void Salvat::printwf(FILE *fout)
{
  int I;
      fprintf(fout, "%d %f %d %f\n", nmax, E_last, k_last, phase_last);
      for(I=0; I<nmax; I++)
        fprintf(fout, "%12.5f %15.8f %15.8f\n", rr[I], P[I], Q[I]);
}
void Salvat::vint(numero *r_in, numero *V_in, int n_in, int atZ)
{
  numero *R_, *RV_;                                   
  int i, nn;                                          
                                                      
  nn=n_in+5;                                          
  R_=new numero[nn];  RV_=new numero[nn];
  R_[0]=0;  RV_[0]=-atZ;                              
  for(i=1; i<nn; i++) if(i<=n_in) {
    R_[i]=r_in[i-1];
    RV_[i]=r_in[i-1]*V_in[i-1];
  }  else  {R_[i]=r_in[n_in-1]+0.001*(i-n_in-1);  RV_[i]=0;}
  vint(R_,RV_,nn);
  delete [] R_;  delete [] RV_;
}
void Salvat::phases(numero *r_in, numero *V_in, int n_in,
                    int atZ,
                    numero Ei, numero Ef, int nE,     
                    int rel, int lmax, numero *ph)    
{                                                     
  numero E;                                           
  int i,l,k;
  vint(r_in, V_in, n_in, atZ);
  for(i=0; i<nE; i++) {
    if(nE==1)  E=Ei;  else  E=Ei+i*(Ef-Ei)/(nE-1);
    if(rel)  for(k=-lmax-1; k<=lmax; k++)  dfree(E,ph[k+lmax+1],k);
    else     for(l=0; l<=lmax; l++)  dfree(E,ph[i*(lmax+1)+l],-l-1);
} }
#endif  
#ifdef jga_Poisson        
#else                     
#define jga_Poisson 1     
                          
                          
numero Poisson(numero r, spline rho)   
{
  if(r==0) r=1e-10;
  return  r*rho.integ(r,rho.b,-1) + rho.integ(1e-10,r);
}
#endif  
#ifdef jga_exchange       
#else                     
#define jga_exchange 1    
                          
                          
numero Vxc_Barth_Hedin_F(numero z)
{
  return (1+z*z*z)*log(1+1/z) + z/2 - z*z -1/3.0;
}
numero Vxc_Barth_Hedin(numero den, numero denup)
{
  if(den<1e-10)  return 0;
  numero rs=1/pow(4*pi/3*den, 1/3.0);
  numero x=denup/den;
  numero alpha0=pow(4/(9*pi), 1/3.0);
  numero a=1/pow(2.0, 1/3.0);
  numero gamma=4*a/(3*(1-a));
  numero cF=0.0254, cP=0.0504;
  numero rF=75,     rP=30;
  numero epsilonFc=-cF*Vxc_Barth_Hedin_F(rs/rF);
  numero epsilonPc=-cP*Vxc_Barth_Hedin_F(rs/rP);
  numero vc=gamma*(epsilonFc-epsilonPc);
  numero epsilon0x=3/(2*pi*alpha0);
  numero epsilonPx=-epsilon0x/rs;
  numero muPx=4*epsilonPx/3;
  numero muPc=-cP*log(1+rP/rs);
  numero AA=muPx+vc;
  numero BB=muPc-vc;
  return  (AA*pow(2*x, 1/3.0) + BB) / 2;  
}
numero Vxc_Slater(numero x)
{
  numero alpha=0.6666666666667;   
  return  -1.5*alpha*pow(3*x/pi,1/3.0);
}
#endif  
#ifdef jga_elements       
#else                     
#define jga_elements 1    
                          
                          
int    element_Z=0;             
numero element_mass=0;          
char   element_symbol[3];       
void element_properties(int Z_)
{
  int    Z=Z_;
  numero mass=0;
  char   symbol[3];
  switch(Z) {
    case  1: strcpy(symbol,"H" );  mass=  1.008;  break;
    case  2: strcpy(symbol,"He");  mass=  4.003;  break;
    case  3: strcpy(symbol,"Li");  mass=  6.941;  break;
    case  4: strcpy(symbol,"Be");  mass=  9.012;  break;
    case  5: strcpy(symbol,"B" );  mass= 10.811;  break;
    case  6: strcpy(symbol,"C" );  mass= 12.011;  break;
    case  7: strcpy(symbol,"N" );  mass= 14.007;  break;
    case  8: strcpy(symbol,"O" );  mass= 15.999;  break;
    case  9: strcpy(symbol,"F" );  mass= 18.998;  break;
    case 10: strcpy(symbol,"Ne");  mass= 20.180;  break;
    case 11: strcpy(symbol,"Na");  mass= 22.990;  break;
    case 12: strcpy(symbol,"Mg");  mass= 24.305;  break;
    case 13: strcpy(symbol,"Al");  mass= 26.982;  break;
    case 14: strcpy(symbol,"Si");  mass= 28.086;  break;
    case 15: strcpy(symbol,"P" );  mass= 30.974;  break;
    case 16: strcpy(symbol,"S" );  mass= 32.066;  break;
    case 17: strcpy(symbol,"Cl");  mass= 35.453;  break;
    case 18: strcpy(symbol,"Ar");  mass= 39.948;  break;
    case 19: strcpy(symbol,"K" );  mass= 39.098;  break;
    case 20: strcpy(symbol,"Ca");  mass= 40.078;  break;
    case 21: strcpy(symbol,"Sc");  mass= 44.956;  break;
    case 22: strcpy(symbol,"Ti");  mass= 47.880;  break;
    case 23: strcpy(symbol,"V" );  mass= 50.942;  break;
    case 24: strcpy(symbol,"Cr");  mass= 51.996;  break;
    case 25: strcpy(symbol,"Mn");  mass= 54.938;  break;
    case 26: strcpy(symbol,"Fe");  mass= 55.847;  break;
    case 27: strcpy(symbol,"Co");  mass= 58.933;  break;
    case 28: strcpy(symbol,"Ni");  mass= 58.693;  break;
    case 29: strcpy(symbol,"Cu");  mass= 63.546;  break;
    case 30: strcpy(symbol,"Zn");  mass= 65.390;  break;
    case 31: strcpy(symbol,"Ga");  mass= 69.723;  break;
    case 32: strcpy(symbol,"Ge");  mass= 72.610;  break;
    case 33: strcpy(symbol,"As");  mass= 74.922;  break;
    case 34: strcpy(symbol,"Se");  mass= 78.960;  break;
    case 35: strcpy(symbol,"Br");  mass= 79.904;  break;
    case 36: strcpy(symbol,"Kr");  mass= 83.800;  break;
    case 37: strcpy(symbol,"Rb");  mass= 85.468;  break;
    case 38: strcpy(symbol,"Sr");  mass= 87.620;  break;
    case 39: strcpy(symbol,"Y" );  mass= 88.906;  break;
    case 40: strcpy(symbol,"Zr");  mass= 91.224;  break;
    case 41: strcpy(symbol,"Nb");  mass= 92.906;  break;
    case 42: strcpy(symbol,"Mo");  mass= 95.940;  break;
    case 43: strcpy(symbol,"Tc");  mass= 97.907;  break;
    case 44: strcpy(symbol,"Ru");  mass=101.070;  break;
    case 45: strcpy(symbol,"Rh");  mass=102.906;  break;
    case 46: strcpy(symbol,"Pd");  mass=106.420;  break;
    case 47: strcpy(symbol,"Ag");  mass=107.868;  break;
    case 48: strcpy(symbol,"Cd");  mass=112.411;  break;
    case 49: strcpy(symbol,"In");  mass=114.818;  break;
    case 50: strcpy(symbol,"Sn");  mass=118.710;  break;
    case 51: strcpy(symbol,"Sb");  mass=121.757;  break;
    case 52: strcpy(symbol,"Te");  mass=127.600;  break;
    case 53: strcpy(symbol,"I" );  mass=126.904;  break;
    case 54: strcpy(symbol,"Xe");  mass=131.290;  break;
    case 55: strcpy(symbol,"Cs");  mass=132.905;  break;
    case 56: strcpy(symbol,"Ba");  mass=137.327;  break;
    case 57: strcpy(symbol,"La");  mass=138.906;  break;
    case 58: strcpy(symbol,"Ce");  mass=140.115;  break;
    case 59: strcpy(symbol,"Pr");  mass=140.908;  break;
    case 60: strcpy(symbol,"Nd");  mass=144.240;  break;
    case 61: strcpy(symbol,"Pm");  mass=144.913;  break;
    case 62: strcpy(symbol,"Sm");  mass=150.360;  break;
    case 63: strcpy(symbol,"Eu");  mass=151.965;  break;
    case 64: strcpy(symbol,"Gd");  mass=157.250;  break;
    case 65: strcpy(symbol,"Tb");  mass=158.925;  break;
    case 66: strcpy(symbol,"Dy");  mass=162.500;  break;
    case 67: strcpy(symbol,"Ho");  mass=164.930;  break;
    case 68: strcpy(symbol,"Er");  mass=167.260;  break;
    case 69: strcpy(symbol,"Tm");  mass=168.934;  break;
    case 70: strcpy(symbol,"Yb");  mass=173.040;  break;
    case 71: strcpy(symbol,"Lu");  mass=174.967;  break;
    case 72: strcpy(symbol,"Hf");  mass=178.490;  break;
    case 73: strcpy(symbol,"Ta");  mass=180.948;  break;
    case 74: strcpy(symbol,"W" );  mass=183.840;  break;
    case 75: strcpy(symbol,"Re");  mass=186.207;  break;
    case 76: strcpy(symbol,"Os");  mass=190.230;  break;
    case 77: strcpy(symbol,"Ir");  mass=192.220;  break;
    case 78: strcpy(symbol,"Pt");  mass=195.080;  break;
    case 79: strcpy(symbol,"Au");  mass=196.967;  break;
    case 80: strcpy(symbol,"Hg");  mass=200.590;  break;
    case 81: strcpy(symbol,"Tl");  mass=204.383;  break;
    case 82: strcpy(symbol,"Pb");  mass=207.200;  break;
    case 83: strcpy(symbol,"Bi");  mass=208.980;  break;
    case 84: strcpy(symbol,"Po");  mass=208.982;  break;
    case 85: strcpy(symbol,"At");  mass=209.987;  break;
    case 86: strcpy(symbol,"Rn");  mass=222.018;  break;
    case 87: strcpy(symbol,"Fr");  mass=223.020;  break;
    case 88: strcpy(symbol,"Ra");  mass=226.025;  break;
    case 89: strcpy(symbol,"Ac");  mass=227.028;  break;
    case 90: strcpy(symbol,"Th");  mass=232.038;  break;
    case 91: strcpy(symbol,"Pa");  mass=231.036;  break;
    case 92: strcpy(symbol,"U" );  mass=238.029;  break;
    case 93: strcpy(symbol,"Np");  mass=237.048;  break;
    case 94: strcpy(symbol,"Pu");  mass=244.064;  break;
    case 95: strcpy(symbol,"Am");  mass=243.061;  break;
    case 96: strcpy(symbol,"Cm");  mass=247.070;  break;
    case 97: strcpy(symbol,"Bk");  mass=247.070;  break;
    case 98: strcpy(symbol,"Cf");  mass=251.080;  break;
    case 99: strcpy(symbol,"Es");  mass=252.083;  break;
    default: strcpy(symbol, "" );  Z=0;           break;
  }
  element_Z=Z;
  element_mass=mass;
  strcpy(element_symbol, symbol);
}
int element_properties(char *symbol)
{
  int Z;
  if(!strcmp(symbol,"H" )) Z= 1; else  if(!strcmp(symbol,"He")) Z= 2; else
  if(!strcmp(symbol,"Li")) Z= 3; else  if(!strcmp(symbol,"Be")) Z= 4; else
  if(!strcmp(symbol,"B" )) Z= 5; else  if(!strcmp(symbol,"C" )) Z= 6; else
  if(!strcmp(symbol,"N" )) Z= 7; else  if(!strcmp(symbol,"O" )) Z= 8; else
  if(!strcmp(symbol,"F" )) Z= 9; else  if(!strcmp(symbol,"Ne")) Z=10; else
  if(!strcmp(symbol,"Na")) Z=11; else  if(!strcmp(symbol,"Mg")) Z=12; else
  if(!strcmp(symbol,"Al")) Z=13; else  if(!strcmp(symbol,"Si")) Z=14; else
  if(!strcmp(symbol,"P" )) Z=15; else  if(!strcmp(symbol,"S" )) Z=16; else
  if(!strcmp(symbol,"Cl")) Z=17; else  if(!strcmp(symbol,"Ar")) Z=18; else
  if(!strcmp(symbol,"K" )) Z=19; else  if(!strcmp(symbol,"Ca")) Z=20; else
  if(!strcmp(symbol,"Sc")) Z=21; else  if(!strcmp(symbol,"Ti")) Z=22; else
  if(!strcmp(symbol,"V" )) Z=23; else  if(!strcmp(symbol,"Cr")) Z=24; else
  if(!strcmp(symbol,"Mn")) Z=25; else  if(!strcmp(symbol,"Fe")) Z=26; else
  if(!strcmp(symbol,"Co")) Z=27; else  if(!strcmp(symbol,"Ni")) Z=28; else
  if(!strcmp(symbol,"Cu")) Z=29; else  if(!strcmp(symbol,"Zn")) Z=30; else
  if(!strcmp(symbol,"Ga")) Z=31; else  if(!strcmp(symbol,"Ge")) Z=32; else
  if(!strcmp(symbol,"As")) Z=33; else  if(!strcmp(symbol,"Se")) Z=34; else
  if(!strcmp(symbol,"Br")) Z=35; else  if(!strcmp(symbol,"Kr")) Z=36; else
  if(!strcmp(symbol,"Rb")) Z=37; else  if(!strcmp(symbol,"Sr")) Z=38; else
  if(!strcmp(symbol,"Y" )) Z=39; else  if(!strcmp(symbol,"Zr")) Z=40; else
  if(!strcmp(symbol,"Nb")) Z=41; else  if(!strcmp(symbol,"Mo")) Z=42; else
  if(!strcmp(symbol,"Tc")) Z=43; else  if(!strcmp(symbol,"Ru")) Z=44; else
  if(!strcmp(symbol,"Rh")) Z=45; else  if(!strcmp(symbol,"Pd")) Z=46; else
  if(!strcmp(symbol,"Ag")) Z=47; else  if(!strcmp(symbol,"Cd")) Z=48; else
  if(!strcmp(symbol,"In")) Z=49; else  if(!strcmp(symbol,"Sn")) Z=50; else
  if(!strcmp(symbol,"Sb")) Z=51; else  if(!strcmp(symbol,"Te")) Z=52; else
  if(!strcmp(symbol,"I" )) Z=53; else  if(!strcmp(symbol,"Xe")) Z=54; else
  if(!strcmp(symbol,"Cs")) Z=55; else  if(!strcmp(symbol,"Ba")) Z=56; else
  if(!strcmp(symbol,"La")) Z=57; else  if(!strcmp(symbol,"Ce")) Z=58; else
  if(!strcmp(symbol,"Pr")) Z=59; else  if(!strcmp(symbol,"Nd")) Z=60; else
  if(!strcmp(symbol,"Pm")) Z=61; else  if(!strcmp(symbol,"Sm")) Z=62; else
  if(!strcmp(symbol,"Eu")) Z=63; else  if(!strcmp(symbol,"Gd")) Z=64; else
  if(!strcmp(symbol,"Tb")) Z=65; else  if(!strcmp(symbol,"Dy")) Z=66; else
  if(!strcmp(symbol,"Ho")) Z=67; else  if(!strcmp(symbol,"Er")) Z=68; else
  if(!strcmp(symbol,"Tm")) Z=69; else  if(!strcmp(symbol,"Yb")) Z=70; else
  if(!strcmp(symbol,"Lu")) Z=71; else  if(!strcmp(symbol,"Hf")) Z=72; else
  if(!strcmp(symbol,"Ta")) Z=73; else  if(!strcmp(symbol,"W" )) Z=74; else
  if(!strcmp(symbol,"Re")) Z=75; else  if(!strcmp(symbol,"Os")) Z=76; else
  if(!strcmp(symbol,"Ir")) Z=77; else  if(!strcmp(symbol,"Pt")) Z=78; else
  if(!strcmp(symbol,"Au")) Z=79; else  if(!strcmp(symbol,"Hg")) Z=80; else
  if(!strcmp(symbol,"Tl")) Z=81; else  if(!strcmp(symbol,"Pb")) Z=82; else
  if(!strcmp(symbol,"Bi")) Z=83; else  if(!strcmp(symbol,"Po")) Z=84; else
  if(!strcmp(symbol,"At")) Z=85; else  if(!strcmp(symbol,"Rn")) Z=86; else
  if(!strcmp(symbol,"Fr")) Z=87; else  if(!strcmp(symbol,"Ra")) Z=88; else
  if(!strcmp(symbol,"Ac")) Z=89; else  if(!strcmp(symbol,"Th")) Z=90; else
  if(!strcmp(symbol,"Pa")) Z=91; else  if(!strcmp(symbol,"U" )) Z=92; else
  if(!strcmp(symbol,"Np")) Z=93; else  if(!strcmp(symbol,"Pu")) Z=94; else
  if(!strcmp(symbol,"Am")) Z=95; else  if(!strcmp(symbol,"Cm")) Z=96; else
  if(!strcmp(symbol,"Bk")) Z=97; else  if(!strcmp(symbol,"Cf")) Z=98; else
  if(!strcmp(symbol,"Es")) Z=99; else                           Z= 0;
  element_properties(Z);
  return Z;
}
void element_orbital(int orb, int &n, int &k)
{
  switch(orb) {
    case  0: k=-1; n=1; break;    case  1: k=-1; n=2; break;
    case  2: k= 1; n=2; break;    case  3: k=-2; n=2; break;
    case  4: k=-1; n=3; break;    case  5: k= 1; n=3; break;
    case  6: k=-2; n=3; break;    case  7: k= 2; n=3; break;
    case  8: k=-3; n=3; break;    case  9: k=-1; n=4; break;
    case 10: k= 1; n=4; break;    case 11: k=-2; n=4; break;
    case 12: k= 2; n=4; break;    case 13: k=-3; n=4; break;
    case 14: k= 3; n=4; break;    case 15: k=-4; n=4; break;
    case 16: k=-1; n=5; break;    case 17: k= 1; n=5; break;
    case 18: k=-2; n=5; break;    case 19: k= 2; n=5; break;
    case 20: k=-3; n=5; break;    case 21: k= 3; n=5; break;
    case 22: k=-4; n=5; break;    case 23: k=-1; n=6; break;
    case 24: k= 1; n=6; break;    case 25: k=-2; n=6; break;
    case 26: k= 2; n=6; break;    case 27: k=-3; n=6; break;
    case 28: k=-1; n=7; break;    case -1: k= 0; n=0; break;
} }
int element_orbital(char *name, int &n, int &k)
{
  int orb=-1, nr=0;       
  if(!strcmp(name,"1s"))     orb= 0;  else   
  if(!strcmp(name,"2s"))     orb= 1;  else   
  if(!strcmp(name,"2p1/2"))  orb= 2;  else   
  if(!strcmp(name,"2p3/2"))  orb= 3;  else   
  if(!strcmp(name,"3s"))     orb= 4;  else   
  if(!strcmp(name,"3p1/2"))  orb= 5;  else   
  if(!strcmp(name,"3p3/2"))  orb= 6;  else   
  if(!strcmp(name,"3d3/2"))  orb= 7;  else   
  if(!strcmp(name,"3d5/2"))  orb= 8;  else   
  if(!strcmp(name,"4s"))     orb= 9;  else   
  if(!strcmp(name,"4p1/2"))  orb=10;  else   
  if(!strcmp(name,"4p3/2"))  orb=11;  else   
  if(!strcmp(name,"4d3/2"))  orb=12;  else   
  if(!strcmp(name,"4d5/2"))  orb=13;  else   
  if(!strcmp(name,"4f5/2"))  orb=14;  else   
  if(!strcmp(name,"4f7/2"))  orb=15;  else   
  if(!strcmp(name,"5s"))     orb=16;  else   
  if(!strcmp(name,"5p1/2"))  orb=17;  else   
  if(!strcmp(name,"5p3/2"))  orb=18;  else   
  if(!strcmp(name,"5d3/2"))  orb=19;  else   
  if(!strcmp(name,"5d5/2"))  orb=20;  else   
  if(!strcmp(name,"5f5/2"))  orb=21;  else   
  if(!strcmp(name,"5f7/2"))  orb=22;  else   
  if(!strcmp(name,"6s"))     orb=23;  else   
  if(!strcmp(name,"6p1/2"))  orb=24;  else   
  if(!strcmp(name,"6p3/2"))  orb=25;  else   
  if(!strcmp(name,"6d3/2"))  orb=26;  else   
  if(!strcmp(name,"6d5/2"))  orb=27;  else   
  if(!strcmp(name,"7s"))     orb=28;  else   
  if(!strcmp(name,"2p"))  {orb= 2;  nr=1;}  else
  if(!strcmp(name,"3p"))  {orb= 5;  nr=1;}  else
  if(!strcmp(name,"3d"))  {orb= 7;  nr=1;}  else
  if(!strcmp(name,"4p"))  {orb=10;  nr=1;}  else
  if(!strcmp(name,"4d"))  {orb=12;  nr=1;}  else
  if(!strcmp(name,"4f"))  {orb=14;  nr=1;}  else
  if(!strcmp(name,"5p"))  {orb=17;  nr=1;}  else
  if(!strcmp(name,"5d"))  {orb=19;  nr=1;}  else
  if(!strcmp(name,"5f"))  {orb=21;  nr=1;}  else
  if(!strcmp(name,"6p"))  {orb=24;  nr=1;}  else
  if(!strcmp(name,"6d"))  {orb=26;  nr=1;}
  element_orbital(orb,n,k);   if(nr) orb=1000;
  return orb;
}
int element_orbital(char *name) {int n,k;  return element_orbital(name,n,k);}
int element_occ(int Z, int orb)
{
  int *a,
  a1[29] ={1,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0},
  a2[29] ={2,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0},
  a3[29] ={2,1,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0},
  a4[29] ={2,2,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0},
  a5[29] ={2,2,1,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0},
  a6[29] ={2,2,2,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0},
  a7[29] ={2,2,2,1,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0},
  a8[29] ={2,2,2,2,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0},
  a9[29] ={2,2,2,3,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0},
  a10[29]={2,2,2,4,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0},
  a11[29]={2,2,2,4,1, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0},
  a12[29]={2,2,2,4,2, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0},
  a13[29]={2,2,2,4,2, 1,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0},
  a14[29]={2,2,2,4,2, 2,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0},
  a15[29]={2,2,2,4,2, 2,1,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0},
  a16[29]={2,2,2,4,2, 2,2,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0},
  a17[29]={2,2,2,4,2, 2,3,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0},
  a18[29]={2,2,2,4,2, 2,4,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0},
  a19[29]={2,2,2,4,2, 2,4,0,0,1, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0},
  a20[29]={2,2,2,4,2, 2,4,0,0,2, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0},
  a21[29]={2,2,2,4,2, 2,4,1,0,2, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0},
  a22[29]={2,2,2,4,2, 2,4,2,0,2, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0},
  a23[29]={2,2,2,4,2, 2,4,3,0,2, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0},
  a24[29]={2,2,2,4,2, 2,4,4,1,1, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0},
  a25[29]={2,2,2,4,2, 2,4,4,1,2, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0},
  a26[29]={2,2,2,4,2, 2,4,4,2,2, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0},
  a27[29]={2,2,2,4,2, 2,4,4,3,2, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0},
  a28[29]={2,2,2,4,2, 2,4,4,4,2, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0},
  a29[29]={2,2,2,4,2, 2,4,4,6,1, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0},
  a30[29]={2,2,2,4,2, 2,4,4,6,2, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0},
  a31[29]={2,2,2,4,2, 2,4,4,6,2, 1,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0},
  a32[29]={2,2,2,4,2, 2,4,4,6,2, 2,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0},
  a33[29]={2,2,2,4,2, 2,4,4,6,2, 2,1,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0},
  a34[29]={2,2,2,4,2, 2,4,4,6,2, 2,2,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0},
  a35[29]={2,2,2,4,2, 2,4,4,6,2, 2,3,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0},
  a36[29]={2,2,2,4,2, 2,4,4,6,2, 2,4,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0},
  a37[29]={2,2,2,4,2, 2,4,4,6,2, 2,4,0,0,0, 0,1,0,0,0, 0,0,0,0,0, 0,0,0,0},
  a38[29]={2,2,2,4,2, 2,4,4,6,2, 2,4,0,0,0, 0,2,0,0,0, 0,0,0,0,0, 0,0,0,0},
  a39[29]={2,2,2,4,2, 2,4,4,6,2, 2,4,1,0,0, 0,2,0,0,0, 0,0,0,0,0, 0,0,0,0},
  a40[29]={2,2,2,4,2, 2,4,4,6,2, 2,4,2,0,0, 0,2,0,0,0, 0,0,0,0,0, 0,0,0,0},
  a41[29]={2,2,2,4,2, 2,4,4,6,2, 2,4,4,0,0, 0,1,0,0,0, 0,0,0,0,0, 0,0,0,0},
  a42[29]={2,2,2,4,2, 2,4,4,6,2, 2,4,4,1,0, 0,1,0,0,0, 0,0,0,0,0, 0,0,0,0},
  a43[29]={2,2,2,4,2, 2,4,4,6,2, 2,4,4,1,0, 0,2,0,0,0, 0,0,0,0,0, 0,0,0,0},
  a44[29]={2,2,2,4,2, 2,4,4,6,2, 2,4,4,3,0, 0,1,0,0,0, 0,0,0,0,0, 0,0,0,0},
  a45[29]={2,2,2,4,2, 2,4,4,6,2, 2,4,4,4,0, 0,1,0,0,0, 0,0,0,0,0, 0,0,0,0},
  a46[29]={2,2,2,4,2, 2,4,4,6,2, 2,4,4,6,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0},
  a47[29]={2,2,2,4,2, 2,4,4,6,2, 2,4,4,6,0, 0,1,0,0,0, 0,0,0,0,0, 0,0,0,0},
  a48[29]={2,2,2,4,2, 2,4,4,6,2, 2,4,4,6,0, 0,2,0,0,0, 0,0,0,0,0, 0,0,0,0},
  a49[29]={2,2,2,4,2, 2,4,4,6,2, 2,4,4,6,0, 0,2,1,0,0, 0,0,0,0,0, 0,0,0,0},
  a50[29]={2,2,2,4,2, 2,4,4,6,2, 2,4,4,6,0, 0,2,2,0,0, 0,0,0,0,0, 0,0,0,0},
  a51[29]={2,2,2,4,2, 2,4,4,6,2, 2,4,4,6,0, 0,2,2,1,0, 0,0,0,0,0, 0,0,0,0},
  a52[29]={2,2,2,4,2, 2,4,4,6,2, 2,4,4,6,0, 0,2,2,2,0, 0,0,0,0,0, 0,0,0,0},
  a53[29]={2,2,2,4,2, 2,4,4,6,2, 2,4,4,6,0, 0,2,2,3,0, 0,0,0,0,0, 0,0,0,0},
  a54[29]={2,2,2,4,2, 2,4,4,6,2, 2,4,4,6,0, 0,2,2,4,0, 0,0,0,0,0, 0,0,0,0},
  a55[29]={2,2,2,4,2, 2,4,4,6,2, 2,4,4,6,0, 0,2,2,4,0, 0,0,0,1,0, 0,0,0,0},
  a56[29]={2,2,2,4,2, 2,4,4,6,2, 2,4,4,6,0, 0,2,2,4,0, 0,0,0,2,0, 0,0,0,0},
  a57[29]={2,2,2,4,2, 2,4,4,6,2, 2,4,4,6,0, 0,2,2,4,1, 0,0,0,2,0, 0,0,0,0},
  a58[29]={2,2,2,4,2, 2,4,4,6,2, 2,4,4,6,2, 0,2,2,4,0, 0,0,0,2,0, 0,0,0,0},
  a59[29]={2,2,2,4,2, 2,4,4,6,2, 2,4,4,6,3, 0,2,2,4,0, 0,0,0,2,0, 0,0,0,0},
  a60[29]={2,2,2,4,2, 2,4,4,6,2, 2,4,4,6,4, 0,2,2,4,0, 0,0,0,2,0, 0,0,0,0},
  a61[29]={2,2,2,4,2, 2,4,4,6,2, 2,4,4,6,5, 0,2,2,4,0, 0,0,0,2,0, 0,0,0,0},
  a62[29]={2,2,2,4,2, 2,4,4,6,2, 2,4,4,6,6, 0,2,2,4,0, 0,0,0,2,0, 0,0,0,0},
  a63[29]={2,2,2,4,2, 2,4,4,6,2, 2,4,4,6,6, 1,2,2,4,0, 0,0,0,2,0, 0,0,0,0},
  a64[29]={2,2,2,4,2, 2,4,4,6,2, 2,4,4,6,6, 1,2,2,4,1, 0,0,0,2,0, 0,0,0,0},
  a65[29]={2,2,2,4,2, 2,4,4,6,2, 2,4,4,6,6, 3,2,2,4,0, 0,0,0,2,0, 0,0,0,0},
  a66[29]={2,2,2,4,2, 2,4,4,6,2, 2,4,4,6,6, 4,2,2,4,0, 0,0,0,2,0, 0,0,0,0},
  a67[29]={2,2,2,4,2, 2,4,4,6,2, 2,4,4,6,6, 5,2,2,4,0, 0,0,0,2,0, 0,0,0,0},
  a68[29]={2,2,2,4,2, 2,4,4,6,2, 2,4,4,6,6, 6,2,2,4,0, 0,0,0,2,0, 0,0,0,0},
  a69[29]={2,2,2,4,2, 2,4,4,6,2, 2,4,4,6,6, 7,2,2,4,0, 0,0,0,2,0, 0,0,0,0},
  a70[29]={2,2,2,4,2, 2,4,4,6,2, 2,4,4,6,6, 8,2,2,4,0, 0,0,0,2,0, 0,0,0,0},
  a71[29]={2,2,2,4,2, 2,4,4,6,2, 2,4,4,6,6, 8,2,2,4,1, 0,0,0,2,0, 0,0,0,0},
  a72[29]={2,2,2,4,2, 2,4,4,6,2, 2,4,4,6,6, 8,2,2,4,2, 0,0,0,2,0, 0,0,0,0},
  a73[29]={2,2,2,4,2, 2,4,4,6,2, 2,4,4,6,6, 8,2,2,4,3, 0,0,0,2,0, 0,0,0,0},
  a74[29]={2,2,2,4,2, 2,4,4,6,2, 2,4,4,6,6, 8,2,2,4,4, 0,0,0,2,0, 0,0,0,0},
  a75[29]={2,2,2,4,2, 2,4,4,6,2, 2,4,4,6,6, 8,2,2,4,4, 1,0,0,2,0, 0,0,0,0},
  a76[29]={2,2,2,4,2, 2,4,4,6,2, 2,4,4,6,6, 8,2,2,4,4, 2,0,0,2,0, 0,0,0,0},
  a77[29]={2,2,2,4,2, 2,4,4,6,2, 2,4,4,6,6, 8,2,2,4,4, 3,0,0,2,0, 0,0,0,0},
  a78[29]={2,2,2,4,2, 2,4,4,6,2, 2,4,4,6,6, 8,2,2,4,4, 5,0,0,1,0, 0,0,0,0},
  a79[29]={2,2,2,4,2, 2,4,4,6,2, 2,4,4,6,6, 8,2,2,4,4, 6,0,0,1,0, 0,0,0,0},
  a80[29]={2,2,2,4,2, 2,4,4,6,2, 2,4,4,6,6, 8,2,2,4,4, 6,0,0,2,0, 0,0,0,0},
  a81[29]={2,2,2,4,2, 2,4,4,6,2, 2,4,4,6,6, 8,2,2,4,4, 6,0,0,2,1, 0,0,0,0},
  a82[29]={2,2,2,4,2, 2,4,4,6,2, 2,4,4,6,6, 8,2,2,4,4, 6,0,0,2,2, 0,0,0,0},
  a83[29]={2,2,2,4,2, 2,4,4,6,2, 2,4,4,6,6, 8,2,2,4,4, 6,0,0,2,2, 1,0,0,0},
  a84[29]={2,2,2,4,2, 2,4,4,6,2, 2,4,4,6,6, 8,2,2,4,4, 6,0,0,2,2, 2,0,0,0},
  a85[29]={2,2,2,4,2, 2,4,4,6,2, 2,4,4,6,6, 8,2,2,4,4, 6,0,0,2,2, 3,0,0,0},
  a86[29]={2,2,2,4,2, 2,4,4,6,2, 2,4,4,6,6, 8,2,2,4,4, 6,0,0,2,2, 4,0,0,0},
  a87[29]={2,2,2,4,2, 2,4,4,6,2, 2,4,4,6,6, 8,2,2,4,4, 6,0,0,2,2, 4,0,0,1},
  a88[29]={2,2,2,4,2, 2,4,4,6,2, 2,4,4,6,6, 8,2,2,4,4, 6,0,0,2,2, 4,0,0,2},
  a89[29]={2,2,2,4,2, 2,4,4,6,2, 2,4,4,6,6, 8,2,2,4,4, 6,0,0,2,2, 4,1,0,2},
  a90[29]={2,2,2,4,2, 2,4,4,6,2, 2,4,4,6,6, 8,2,2,4,4, 6,0,0,2,2, 4,2,0,2},
  a91[29]={2,2,2,4,2, 2,4,4,6,2, 2,4,4,6,6, 8,2,2,4,4, 6,2,0,2,2, 4,1,0,2},
  a92[29]={2,2,2,4,2, 2,4,4,6,2, 2,4,4,6,6, 8,2,2,4,4, 6,3,0,2,2, 4,1,0,2},
  a93[29]={2,2,2,4,2, 2,4,4,6,2, 2,4,4,6,6, 8,2,2,4,4, 6,4,0,2,2, 4,1,0,2},
  a94[29]={2,2,2,4,2, 2,4,4,6,2, 2,4,4,6,6, 8,2,2,4,4, 6,6,0,2,2, 4,0,0,2},
  a95[29]={2,2,2,4,2, 2,4,4,6,2, 2,4,4,6,6, 8,2,2,4,4, 6,6,1,2,2, 4,0,0,2};
  if(Z== 1) a=a1;  else if(Z== 2) a=a2;  else if(Z== 3) a=a3;  else
  if(Z== 4) a=a4;  else if(Z== 5) a=a5;  else if(Z== 6) a=a6;  else
  if(Z== 7) a=a7;  else if(Z== 8) a=a8;  else if(Z== 9) a=a9;  else
  if(Z==10) a=a10; else
  if(Z==11) a=a11; else if(Z==12) a=a12; else if(Z==13) a=a13; else
  if(Z==14) a=a14; else if(Z==15) a=a15; else if(Z==16) a=a16; else
  if(Z==17) a=a17; else if(Z==18) a=a18; else if(Z==19) a=a19; else
  if(Z==20) a=a20; else
  if(Z==21) a=a21; else if(Z==22) a=a22; else if(Z==23) a=a23; else
  if(Z==24) a=a24; else if(Z==25) a=a25; else if(Z==26) a=a26; else
  if(Z==27) a=a27; else if(Z==28) a=a28; else if(Z==29) a=a29; else
  if(Z==30) a=a30; else
  if(Z==31) a=a31; else if(Z==32) a=a32; else if(Z==33) a=a33; else
  if(Z==34) a=a34; else if(Z==35) a=a35; else if(Z==36) a=a36; else
  if(Z==37) a=a37; else if(Z==38) a=a38; else if(Z==39) a=a39; else
  if(Z==40) a=a40; else
  if(Z==41) a=a41; else if(Z==42) a=a42; else if(Z==43) a=a43; else
  if(Z==44) a=a44; else if(Z==45) a=a45; else if(Z==46) a=a46; else
  if(Z==47) a=a47; else if(Z==48) a=a48; else if(Z==49) a=a49; else
  if(Z==50) a=a50; else
  if(Z==51) a=a51; else if(Z==52) a=a52; else if(Z==53) a=a53; else
  if(Z==54) a=a54; else if(Z==55) a=a55; else if(Z==56) a=a56; else
  if(Z==57) a=a57; else if(Z==58) a=a58; else if(Z==59) a=a59; else
  if(Z==60) a=a60; else
  if(Z==61) a=a61; else if(Z==62) a=a62; else if(Z==63) a=a63; else
  if(Z==64) a=a64; else if(Z==65) a=a65; else if(Z==66) a=a66; else
  if(Z==67) a=a67; else if(Z==68) a=a68; else if(Z==69) a=a69; else
  if(Z==70) a=a70; else
  if(Z==71) a=a71; else if(Z==72) a=a72; else if(Z==73) a=a73; else
  if(Z==74) a=a74; else if(Z==75) a=a75; else if(Z==76) a=a76; else
  if(Z==77) a=a77; else if(Z==78) a=a78; else if(Z==79) a=a79; else
  if(Z==80) a=a80; else
  if(Z==81) a=a81; else if(Z==82) a=a82; else if(Z==83) a=a83; else
  if(Z==84) a=a84; else if(Z==85) a=a85; else if(Z==86) a=a86; else
  if(Z==87) a=a87; else if(Z==88) a=a88; else if(Z==89) a=a89; else
  if(Z==90) a=a90; else
  if(Z==91) a=a91; else if(Z==92) a=a92; else if(Z==93) a=a93; else
  if(Z==94) a=a94; else if(Z==95) a=a95;
  return a[orb];
}
#endif  
#ifdef jga_Rho            
#else                     
#define jga_Rho 1         
                          
#define Rho_Zmax  95           
#define Rho_imax 231           
numero Rho_radius[231]={  
	              0,  0.000158461291,  0.000166585785,  0.000175126814, 
	 0.000184105753,  0.000193545056,  0.000203468328,  0.000213900377, 
	 0.000224867283,  0.000236396474,  0.000248516793,  0.000261258509, 
	 0.000274653517,  0.000288735289,  0.000303539069,  0.000319101848, 
	 0.000335462566,  0.000352662086,  0.000370743481,  0.000389751891, 
	 0.000409734901,  0.000430742453,   0.00045282708,  0.000476044021, 
	 0.000500451308,  0.000526110001,  0.000553084246,  0.000581441505, 
	 0.000611252617,  0.000642592204,  0.000675538613,  0.000710174267, 
	 0.000746585662,  0.000784863951,  0.000825104769,  0.000867408817, 
	 0.000911881798,  0.000958634948,   0.00100778521,   0.00105945545, 
	  0.00111377495,   0.00117087935,    0.0012309117,   0.00129402184, 
	  0.00136036775,   0.00143011531,   0.00150343892,   0.00158052193, 
	  0.00166155689,   0.00174674683,   0.00183630444,   0.00193045381, 
	  0.00202943012,    0.0021334812,   0.00224286737,   0.00235786149, 
	  0.00247875159,    0.0026058401,   0.00273944438,   0.00287989853, 
	  0.00302755414,   0.00318278023,   0.00334596471,   0.00351751619, 
	  0.00369786308,   0.00388745661,   0.00408677058,   0.00429630373, 
	  0.00451658014,   0.00474815024,   0.00499159284,   0.00524751749, 
	  0.00551656308,   0.00579940341,   0.00609674538,   0.00640933216, 
	   0.0067379456,    0.0070834076,    0.0074465815,   0.00782837626, 
	  0.00822974555,   0.00865169335,   0.00909527484,   0.00956160016, 
	   0.0100518344,    0.0105672022,    0.0111089936,    0.0116785653, 
	   0.0122773368,    0.0129068103,    0.0135685559,    0.0142642306, 
	    0.014995574,    0.0157644134,    0.0165726729,    0.0174223706, 
	   0.0183156356,    0.0192546975,    0.0202419069,    0.0212797318, 
	   0.0223707668,    0.0235177409,    0.0247235224,    0.0259911232, 
	   0.0273237173,    0.0287246332,    0.0301973764,    0.0317456312, 
	   0.0333732627,    0.0350843482,    0.0368831605,    0.0387741998, 
	   0.0407621972,    0.0428521186,    0.0450491942,    0.0473589152, 
	   0.0497870594,    0.0523396954,    0.0550232083,    0.0578443073, 
	   0.0608100519,    0.0639278516,    0.0672054961,    0.0706511959, 
	   0.0742735639,    0.0780816525,    0.0820849836,    0.0862935707, 
	   0.0907179341,    0.0953691453,      0.10025882,     0.105399199, 
	    0.110803142,     0.116484128,     0.122456402,     0.128734887, 
	    0.135335267,     0.142274037,     0.149568588,     0.157237142, 
	    0.165298864,     0.173773915,     0.182683483,     0.192049876, 
	    0.201896474,     0.212247923,     0.223130122,     0.234570235, 
	    0.246596918,      0.25924021,     0.272531748,     0.286504745, 
	    0.301194161,     0.316636711,      0.33287102,     0.349937677, 
	    0.367879361,     0.386740953,      0.40656957,     0.427414864, 
	     0.44932887,     0.472366452,      0.49658522,     0.522045672, 
	    0.548811555,     0.576949716,     0.606530547,     0.637628019, 
	    0.670319915,     0.704687953,     0.740818083,     0.778800607, 
	    0.818730593,     0.860707819,     0.904837251,     0.951229215, 
	    0.999999821,      1.05127084,      1.10517073,        1.161834, 
	     1.22140253,      1.28402519,      1.34985852,      1.41906726, 
	     1.49182439,      1.56831193,      1.64872086,      1.73325264, 
	      1.8221184,      1.91554046,      2.01375222,      2.11699963, 
	      2.2255404,      2.33964634,      2.45960259,       2.5857091, 
	     2.71828127,      2.85765052,      3.00416541,       3.1581924, 
	     3.32011628,      3.49034238,      3.66929579,      3.85742474, 
	     4.05519915,       4.2631135,      4.48168802,      4.71146917, 
	     4.95303154,       5.2069788,      5.47394609,      5.75460148, 
	     6.04964638,      6.35981846,      6.68589306,      7.02868605, 
	     7.38905478,      7.76789951,      8.16616821,      8.58485699, 
	     9.02501202,      9.48773384,      9.97418022,       10.485568, 
	     11.0231743,      11.5883436,      12.1824923,      12.8071012, 
	     13.4637346,      14.1540356,      14.8797293};
numero Rho_density_data[21945]={  
/* Z= 1 -> H */
	              0,  1.00456127e-07,  1.11019091e-07,  1.22692654e-07, 
	 1.35593552e-07,   1.4985082e-07,  1.65607034e-07,  1.83019779e-07, 
	 2.02263166e-07,  2.23529611e-07,  2.47031778e-07,  2.73004673e-07, 
	 3.01707928e-07,  3.33428545e-07,  3.68483626e-07,  4.07223638e-07, 
	 4.50035799e-07,  4.97348026e-07,  5.49633285e-07,  6.07413995e-07, 
	 6.71267628e-07,  7.41832366e-07,  8.19813124e-07,  9.05989168e-07, 
	 1.00122134e-06,  1.10646101e-06,  1.22275935e-06,  1.35127766e-06, 
	 1.49329981e-06,  1.65024369e-06,  1.82367626e-06,  2.01532907e-06, 
	 2.22711492e-06,  2.46114769e-06,  2.71976273e-06,  3.00554029e-06, 
	 3.32133141e-06,  3.67028611e-06,  4.05588389e-06,  4.48196988e-06, 
	 4.95279119e-06,  5.47304126e-06,  6.04790375e-06,  6.68310622e-06, 
	 7.38497465e-06,  8.16049942e-06,  9.01739986e-06,  9.96420476e-06, 
	 1.10103356e-05,  1.21661978e-05,  1.34432839e-05,  1.48542904e-05, 
	 1.64132362e-05,  1.81356081e-05,  2.00385093e-05,  2.21408263e-05, 
	 2.44634157e-05,  2.70293131e-05,  2.98639497e-05,  3.29954128e-05, 
	  3.6454705e-05,   4.0276067e-05,  4.44972939e-05,  4.91601131e-05, 
	 5.43105925e-05,  5.99995692e-05,  6.62831735e-05,   7.3223353e-05, 
	 8.08884579e-05,    8.935394e-05,  9.87030435e-05,  0.000109027627, 
	 0.000120429017,  0.000133019028,  0.000146920967,  0.000162270881, 
	 0.000179218754,  0.000197930014,  0.000218587084,  0.000241391055, 
	 0.000266563613,   0.00029434904,  0.000325016677,  0.000358863181, 
	 0.000396215415,  0.000437433482,  0.000482913922,  0.000533093407, 
	 0.000588452618,  0.000649520778,  0.000716880138,  0.000791171391, 
	 0.000873099256,  0.000963438768,   0.00106304185,   0.00117284502, 
	  0.00129387714,   0.00142726838,   0.00157425995,   0.00173621427, 
	  0.00191462599,   0.00211113482,   0.00232753903,   0.00256580836, 
	  0.00282810186,   0.00311678252,   0.00343443733,   0.00378389517, 
	  0.00416824874,   0.00459087733,   0.00505547039,   0.00556605263, 
	  0.00612701382,    0.0067431354,    0.0074196225,   0.00816213712, 
	  0.00897683296,   0.00987038855,    0.0108500477,    0.0119236596, 
	   0.0130997151,    0.0143873906,    0.0157965925,    0.0173379928, 
	   0.0190230776,    0.0208641831,    0.0228745434,    0.0250683185, 
	   0.0274606328,    0.0300676059,    0.0329063721,    0.0359951034, 
	   0.0393530056,    0.0430003256,    0.0469583236,     0.051249247, 
	   0.0558962598,    0.0609233938,    0.0663554221,    0.0722177401, 
	   0.0785361975,    0.0853368938,    0.0926459357,      0.10048914, 
	    0.108891688,     0.117877759,     0.127470016,     0.137689099, 
	    0.148553059,     0.160076678,      0.17227073,     0.185141161, 
	    0.198688269,     0.212905705,     0.227779493,     0.243287057, 
	    0.259396046,     0.276063293,     0.293233812,     0.310839623, 
	     0.32879892,     0.347015172,     0.365376502,     0.383755177, 
	    0.402007371,     0.419973403,     0.437478274,      0.45433253, 
	    0.470333964,     0.485269397,     0.498917788,     0.511053145, 
	    0.521448731,     0.529881775,     0.536138415,     0.540019751, 
	    0.541347742,     0.539971888,     0.535775244,     0.528681159, 
	    0.518658578,      0.50572741,      0.48996228,     0.471494913, 
	    0.450515002,     0.427268982,     0.402056634,     0.375225365, 
	    0.347162366,     0.318284601,     0.289026827,     0.259828627, 
	    0.231120169,     0.203308254,     0.176762715,     0.151804343, 
	    0.128695041,     0.107630543,    0.0887365714,    0.0720684081, 
	   0.0576139838,    0.0453001708,    0.0350018516,    0.0265529528, 
	   0.0197586175,    0.0144076366,    0.0102841826,   0.00717815617, 
	  0.00489353901,   0.00325443409,    0.0021087185,    0.0013294518, 
	  0.00081438123,  0.000483997515,  0.000278641004,  0.000155140442, 
	 8.33947197e-05,  4.32019078e-05,  2.15275886e-05,  1.02979602e-05, 
	 4.71912699e-06,  2.06715026e-06,  8.63534012e-07,  3.43186997e-07, 
	 1.29425828e-07,   4.6194522e-08,  1.55605235e-08,  4.93230257e-09, 
	 1.46668255e-09,  4.07843564e-10,  1.05699609e-10, 
/* Z= 2 -> He */
	              0,  1.13589454e-06,  1.25530335e-06,  1.38726273e-06, 
	 1.53309099e-06,  1.69424561e-06,  1.87233672e-06,  2.06914365e-06, 
	 2.28663271e-06,  2.52697646e-06,  2.79257597e-06,  3.08608355e-06, 
	 3.41043074e-06,  3.76885646e-06,  4.16493958e-06,  4.60263436e-06, 
	 5.08631092e-06,  5.62079595e-06,  6.21142499e-06,  6.86409112e-06, 
	 7.58530723e-06,  8.38226697e-06,  9.26292068e-06,  1.02360509e-05, 
	  1.1311362e-05,  1.24995722e-05,   1.3812526e-05,  1.52633074e-05, 
	 1.68663737e-05,   1.8637691e-05,  2.05949018e-05,  2.27574892e-05, 
	 2.51469864e-05,  2.77871677e-05,   3.0704301e-05,  3.39273975e-05, 
	 3.74885058e-05,  4.14230199e-05,  4.57700335e-05,  5.05727185e-05, 
	 5.58787651e-05,  6.17408296e-05,  6.82170576e-05,  7.53716813e-05, 
	 8.32756123e-05,  9.20071398e-05,  0.000101652717,  0.000112307811, 
	 0.000124077793,  0.000137079012,  0.000151439875,  0.000167302176, 
	 0.000184822376,  0.000204173179,  0.000225545198,  0.000249148754, 
	 0.000275215978,  0.000304002984,  0.000335792254,  0.000370895577, 
	 0.000409656757,  0.000452455075,  0.000499708811,  0.000551879173, 
	 0.000609474839,  0.000673056464,  0.000743242214,  0.000820713409, 
	 0.000906220812,    0.0010005919,    0.0011047381,   0.00121966354, 
	  0.00134647405,   0.00148638745,   0.00164074451,   0.00181102089, 
	  0.00199884083,   0.00220599119,    0.0024344374,   0.00268634153, 
	  0.00296408031,   0.00327026588,   0.00360776903,   0.00397974299, 
	  0.00438965019,   0.00484129321,   0.00533884251,   0.00588687463, 
	  0.00649040798,   0.00715494249,   0.00788650475,   0.00869169179, 
	  0.00957772788,    0.0105525143,    0.0116246901,    0.0128036942, 
	   0.0140998363,    0.0155243669,    0.0170895569,    0.0188087728, 
	   0.0206965767,    0.0227688085,    0.0250426866,    0.0275369156, 
	   0.0302717797,    0.0332692713,    0.0365532003,    0.0401492976, 
	   0.0440853685,    0.0483913869,    0.0530996323,    0.0582448132, 
	    0.063864179,    0.0699976459,    0.0766878873,    0.0839804485, 
	   0.0919238254,     0.100569502,     0.109972037,     0.120189048, 
	    0.131281212,     0.143312156,     0.156348407,     0.170459226, 
	    0.185716346,     0.202193692,     0.219966978,     0.239113286, 
	    0.259710371,     0.281836092,     0.305567503,     0.330979884, 
	    0.358145654,     0.387133032,     0.418004692,     0.450816005, 
	    0.485613376,     0.522432148,     0.561294496,     0.602207065, 
	    0.645158589,      0.69011718,     0.737027764,       0.7858091, 
	    0.836351335,     0.888513088,     0.942119002,     0.996957183, 
	     1.05277741,      1.10928953,      1.16616201,      1.22302175, 
	     1.27945435,      1.33500564,      1.38918328,      1.44146061, 
	     1.49128115,       1.5380646,      1.58121407,      1.62012458, 
	     1.65419304,      1.68282926,      1.70546842,      1.72158372, 
	     1.73070037,      1.73240936,      1.72638035,      1.71237636, 
	     1.69026423,      1.66002572,      1.62176561,      1.57571781, 
	     1.52224755,        1.461851,       1.3951503,      1.32288587, 
	     1.24590361,      1.16513884,      1.08159697,     0.996330559, 
	    0.910414219,     0.824918032,     0.740880013,     0.659279406, 
	    0.581010818,     0.506861448,     0.437491953,     0.373421729, 
	    0.315019578,     0.262500376,     0.215927616,     0.175222009, 
	    0.140175432,     0.110469371,    0.0856969506,    0.0653869659, 
	    0.049028568,     0.036095228,    0.0260665063,    0.0184466112, 
	   0.0127790077,   0.00865659025,   0.00572748017,   0.00369674386, 
	  0.00232465775,   0.00142232596,  0.000845524773,  0.000487636047, 
	 0.000272414298,  0.000147169005,  7.67549645e-05,  3.85756684e-05, 
	 1.86470625e-05,  8.65224047e-06,  3.84551959e-06,  1.63354457e-06, 
	 6.61682236e-07,  2.54949157e-07,   9.3203937e-08,  3.22425073e-08, 
	 1.05249232e-08,  3.23243388e-09,  9.31168087e-10,  2.50797549e-10, 
	  6.2945725e-11,  1.46707004e-11,  3.16385134e-12,  6.28997017e-13, 
	 1.14838464e-13,   1.9179157e-14,  2.91834169e-15, 
/* Z= 3 -> Li */
	              0,  4.37118069e-06,   4.8305501e-06,  5.33818229e-06, 
	 5.89914453e-06,  6.51903747e-06,  7.20404933e-06,  7.96101631e-06, 
	 8.79749405e-06,  9.72182897e-06,  1.07432443e-05,   1.1871929e-05, 
	 1.31191418e-05,  1.44973219e-05,  1.60202126e-05,  1.77029942e-05, 
	 1.95624452e-05,  2.16170974e-05,  2.38874236e-05,  2.63960428e-05, 
	 2.91679426e-05,  3.22307278e-05,  3.56148885e-05,  3.93541159e-05, 
	 4.34856156e-05,  4.80504859e-05,  5.30941361e-05,  5.86667084e-05, 
	 6.48235946e-05,  7.16259674e-05,  7.91413913e-05,  8.74444959e-05, 
	 9.66176958e-05,  0.000106751992,  0.000117947893,  0.000130316374, 
	 0.000143980011,  0.000159074087,   0.00017574802,  0.000194166743, 
	 0.000214512402,  0.000236986001,  0.000261809473,  0.000289227755, 
	  0.00031951128,  0.000352958421,  0.000389898516,  0.000430695014, 
	 0.000475748908,  0.000525502721,  0.000580444525,  0.000641112914, 
	 0.000708101958,  0.000782066782,  0.000863729976,  0.000953888288, 
	  0.00105342036,   0.00116329477,   0.00128457928,    0.0014184505, 
	  0.00156620552,   0.00172927312,   0.00190922769,   0.00210780324, 
	  0.00232690992,   0.00256865029,   0.00283533987,   0.00312952697, 
	  0.00345401559,   0.00381189166,   0.00420654891,   0.00464171916, 
	  0.00512150582,   0.00565041834,   0.00623341277,   0.00687593129, 
	  0.00758395111,   0.00836403295,   0.00922337733,    0.0101698805, 
	   0.0112122009,    0.0123598287,    0.0136231612,    0.0150135821, 
	   0.0165435486,    0.0182266925,    0.0200779103,    0.0221134815, 
	   0.0243511815,    0.0268104076,    0.0295123067,    0.0324799232, 
	   0.0357383527,    0.0393148847,      0.04323918,    0.0475434512, 
	   0.0522626415,    0.0574346073,    0.0631003156,    0.0693040714, 
	   0.0760936886,    0.0835207254,    0.0916406661,      0.10051316, 
	    0.110202193,     0.120776303,     0.132308736,     0.144877598, 
	    0.158566013,     0.173462197,     0.189659536,     0.207256526, 
	    0.226356849,     0.247069106,     0.269506723,     0.293787479, 
	    0.320033342,     0.348369628,     0.378924429,     0.411827773, 
	    0.447210491,     0.485202938,     0.525933564,      0.56952697, 
	    0.616102219,     0.665770113,      0.71863091,     0.774771333, 
	    0.834261179,     0.897149861,     0.963462591,      1.03319609, 
	     1.10631442,      1.18274379,      1.26236844,      1.34502542, 
	     1.43049955,      1.51851881,      1.60875022,      1.70079494, 
	     1.79418564,      1.88838315,      1.98277557,      2.07667708, 
	     2.16932917,      2.25990391,      2.34750891,      2.43119287, 
	     2.50995684,      2.58276439,      2.64855766,      2.70627236, 
	     2.75485945,       2.7933054,       2.8206563,      2.83604407, 
	     2.83871174,      2.82804012,       2.8035748,      2.76504874, 
	     2.71240592,      2.64581919,      2.56570244,      2.47272015, 
	     2.36778736,      2.25206351,      2.12693787,      1.99400628, 
	     1.85504138,      1.71195269,      1.56674206,      1.42145276, 
	     1.27811456,      1.13868833,      1.00501025,     0.878740668, 
	    0.761316955,     0.653916001,     0.557426214,     0.472430438, 
	    0.399202317,     0.337713957,     0.287655681,     0.248465836, 
	    0.219368368,     0.199415877,     0.187534973,     0.182571426, 
	    0.183332443,     0.188624233,     0.197283581,     0.208203018, 
	    0.220349818,     0.232779726,     0.244646654,     0.255209774, 
	    0.263839036,     0.270020127,     0.273359478,     0.273588419, 
	    0.270566821,     0.264284313,     0.254858524,     0.242528975, 
	    0.227645829,     0.210653082,     0.192066744,     0.172448367, 
	    0.152375966,     0.132413641,     0.113082796,    0.0948366374, 
	   0.0780402049,    0.0629576147,    0.0497471541,    0.0384642333, 
	   0.0290716179,    0.0214552432,    0.0154437311,    0.0108292922, 
	  0.00738792494,   0.00489705475,   0.00314934971,    0.0019621423, 
	  0.00118244125,  0.000688092725,  0.000385987078,  0.000208333237, 
	 0.000107984568,  5.36411208e-05,  2.54821152e-05, 
/* Z= 4 -> Be */
	              0,  1.12235603e-05,  1.24026174e-05,  1.37054913e-05, 
	 1.51451768e-05,  1.67360322e-05,  1.84939199e-05,   2.0436366e-05, 
	 2.25827334e-05,  2.49544155e-05,  2.75750463e-05,  3.04707355e-05, 
	  3.3670327e-05,  3.72056893e-05,  4.11120272e-05,  4.54282235e-05, 
	 5.01972427e-05,  5.54665385e-05,  6.12885269e-05,  6.77211137e-05, 
	 7.48282546e-05,  8.26805917e-05,  9.13561526e-05,  0.000100941121, 
	 0.000111530659,  0.000123229896,  0.000136154908,  0.000150433916, 
	 0.000166208454,  0.000183634897,  0.000202885814,  0.000224151852, 
	 0.000247643417,  0.000273592857,  0.000302256696,  0.000333918084, 
	 0.000368889596,  0.000407516258,  0.000450178923,  0.000497297908, 
	  0.00054933707,  0.000606808288,  0.000670276408,   0.00074036466, 
	 0.000817760592,  0.000903222768,  0.000997587922,   0.00110177882, 
	  0.00121681322,   0.00134381349,   0.00148401712,   0.00163878838, 
	  0.00180963136,   0.00199820357,   0.00220633158,   0.00243602879, 
	  0.00268951268,    0.0029692261,   0.00327785918,   0.00361837493, 
	   0.0039940346,   0.00440842891,   0.00486550853,   0.00536962086, 
	   0.0059255478,   0.00653854758,   0.00721440231,   0.00795946736, 
	  0.00878072437,   0.00968584698,    0.0106832562,     0.011782201, 
	   0.0129928291,    0.0143262716,    0.0157947373,    0.0174116064, 
	   0.0191915352,    0.0211505778,    0.0233063046,    0.0256779343, 
	   0.0282864869,    0.0311549231,    0.0343083218,    0.0377740525, 
	   0.0415819585,    0.0457645655,    0.0503572971,    0.0553986803, 
	    0.060930606,    0.0669985637,    0.0736519024,    0.0809441209, 
	   0.0889331102,    0.0976814851,     0.107256837,     0.117732093, 
	    0.129185721,     0.141702116,     0.155371845,      0.17029193, 
	    0.186566085,     0.204304948,      0.22362633,     0.244655281, 
	    0.267524213,     0.292372912,     0.319348484,     0.348605156, 
	    0.380304039,     0.414612651,     0.451704293,     0.491757363, 
	     0.53495425,     0.581480086,     0.631521523,     0.685264409, 
	    0.742892087,     0.804582715,     0.870506346,     0.940821767, 
	     1.01567256,      1.09518301,      1.17945325,      1.26855481, 
	     1.36252367,      1.46135569,       1.5649991,      1.67334831, 
	     1.78623664,      1.90342879,      2.02461457,      2.14940095, 
	      2.2773056,      2.40775061,      2.54005814,      2.67344475, 
	     2.80702019,      2.93978572,      3.07063532,      3.19836044, 
	     3.32165647,      3.43913317,      3.54932737,       3.6507225, 
	     3.74176788,      3.82090712,      3.88660479,      3.93738174, 
	     3.97185063,      3.98875499,      3.98701024,      3.96574521, 
	     3.92434144,      3.86247325,      3.78014088,      3.67770028, 
	     3.55588436,      3.41581464,      3.25900388,      3.08734488, 
	     2.90308762,      2.70880103,      2.50732374,      2.30170059, 
	     2.09510756,      1.89076924,      1.69186831,      1.50145376, 
	     1.32234931,      1.15706742,      1.00773227,     0.876016557, 
	    0.763094068,     0.669611156,     0.595678687,     0.540884018, 
	     0.50432235,     0.484645575,     0.480124861,     0.488724113, 
	    0.508179426,     0.536081493,     0.569956541,     0.607343674, 
	    0.645865917,     0.683293641,     0.717599988,      0.74700731, 
	    0.770025611,     0.785481334,      0.79253757,     0.790704727, 
	    0.779841065,     0.760142148,     0.732119501,     0.696568251, 
	    0.654524326,     0.607212126,     0.555986106,      0.50226754, 
	    0.447480291,      0.39298895,     0.340042859,     0.289728671, 
	    0.242934436,     0.200326711,      0.16234149,     0.129189074, 
	     0.10087131,    0.0772093907,    0.0578791834,    0.0424509794, 
	   0.0304300357,    0.0212948974,    0.0145307416,   0.00965585746, 
	  0.00624031154,   0.00391673716,     0.002383983,   0.00140496797, 
	 0.000800392008,  0.000440010946,   0.00023300288,  0.000118622753, 
	 5.79445805e-05,  2.71006902e-05,  1.21090015e-05,  5.15686543e-06, 
	 2.08809661e-06,  8.01840258e-07,  2.91224268e-07, 
/* Z= 5 -> B */
	              0,  2.28794652e-05,  2.52819827e-05,   2.7936665e-05, 
	 3.08699637e-05,  3.41110936e-05,  3.76923381e-05,  4.16493567e-05, 
	 4.60215488e-05,  5.08524281e-05,  5.61900706e-05,  6.20875871e-05, 
	 6.86036365e-05,  7.58030219e-05,   8.3757317e-05,  9.25455897e-05, 
	 0.000102255151,  0.000112982467,  0.000124834056,  0.000137927564, 
	 0.000152392924,  0.000168373648,  0.000186028192,  0.000205531542, 
	 0.000227076933,  0.000250877754,  0.000277169573,  0.000306212518, 
	 0.000338293787,  0.000373730436,  0.000412872498,  0.000456106296, 
	 0.000503858435,  0.000556599523,  0.000614849268,  0.000679180957, 
	 0.000750227424,  0.000828686869,  0.000915329787,   0.00101100607, 
	  0.00111665356,   0.00123330671,   0.00136210618,   0.00150431064, 
	  0.00166130741,   0.00183462666,   0.00202595559,   0.00223715371, 
	  0.00247027073,   0.00272756605,   0.00301152887,   0.00332490192, 
	  0.00367070711,   0.00405227207,     0.004473262,   0.00493771117, 
	  0.00545006292,   0.00601520482,   0.00663851714,   0.00732591748, 
	  0.00808391441,   0.00891966559,   0.00984103885,    0.0108566787, 
	   0.0119760865,     0.013209695,    0.0145689566,    0.0160664488, 
	   0.0177159607,    0.0195326228,    0.0215330217,    0.0237353314, 
	   0.0261594597,    0.0288272016,    0.0317624025,    0.0349911526, 
	   0.0385419503,     0.042445939,    0.0467371196,    0.0514525808, 
	    0.056632746,    0.0623216666,    0.0685672686,    0.0754216984, 
	   0.0829415843,    0.0911884233,     0.100228868,     0.110135138, 
	    0.120985359,     0.132863939,     0.145861953,     0.160077572, 
	    0.175616369,     0.192591786,     0.211125404,     0.231347367, 
	     0.25339669,     0.277421415,     0.303579003,     0.332036346, 
	    0.362969905,     0.396565646,     0.433018923,      0.47253409, 
	    0.515324056,     0.561609745,     0.611618996,     0.665585399, 
	    0.723747134,     0.786344647,     0.853619039,     0.925809085, 
	     1.00314832,      1.08586144,      1.17416024,      1.26823938, 
	      1.3682704,      1.47439635,      1.58672583,      1.70532441, 
	      1.8302089,      1.96133745,      2.09860229,      2.24181914, 
	     2.39071918,       2.5449388,      2.70401025,      2.86735344, 
	       3.034266,      3.20391679,      3.37534022,      3.54743075, 
	     3.71894145,      3.88848472,      4.05453491,      4.21543789, 
	     4.36942101,      4.51460981,      4.64905167,       4.7707386, 
	     4.87764263,      4.96775293,      5.03911877,      5.08989477, 
	     5.11839628,      5.12314987,      5.10295343,      5.05692577, 
	     4.98456669,      4.88580036,      4.76102114,      4.61112309, 
	     4.43752527,      4.24217367,      4.02753496,      3.79656816, 
	     3.55267763,        3.299649,      3.04156494,      2.78270507, 
	      2.5274322,      2.28006744,      2.04476094,      1.82536161, 
	     1.62529337,      1.44744205,      1.29405975,      1.16669226, 
	     1.06613076,     0.992394924,     0.944742262,      0.92170918, 
	    0.921177208,     0.940462112,     0.976421416,      1.02557361, 
	      1.0842247,      1.14859617,      1.21494865,      1.27969861, 
	     1.33952415,      1.39145756,      1.43296206,       1.4619931, 
	     1.47704196,      1.47716093,      1.46196938,       1.4316411, 
	     1.38687205,      1.32883096,      1.25909138,      1.17955256, 
	     1.09234738,     0.999746323,     0.904058278,     0.807535708, 
	    0.712287724,     0.620205402,     0.532902896,     0.451675922, 
	    0.377480209,     0.310928017,     0.252302974,     0.201589778, 
	    0.158516213,     0.122603163,    0.0932189971,    0.0696340278, 
	   0.0510715432,    0.0367525741,    0.0259323418,     0.017927289, 
	   0.0121326819,   0.00803139713,   0.00519527216,   0.00328066759, 
	  0.00202006195,   0.00121138082,  0.000706519233,  0.000400178134, 
	 0.000219770736,  0.000116819407,  5.99887753e-05,  2.96999751e-05, 
	 1.41460578e-05,  6.46714489e-06,  2.83097779e-06,  1.18359128e-06, 
	 4.71357936e-07,  1.78309506e-07,  6.38870645e-08, 
/* Z= 6 -> C */
	              0,   4.0714538e-05,  4.49878207e-05,  4.97093715e-05, 
	  5.4926164e-05,  6.06901049e-05,  6.70585214e-05,  7.40947435e-05, 
	 8.18687331e-05,  9.04577537e-05,  9.99471522e-05,   0.00011043121, 
	 0.000122014048,  0.000134810674,  0.000148948093,  0.000164566576, 
	 0.000181821059,  0.000200882627,  0.000221940179,  0.000245202391, 
	 0.000270899589,  0.000299286214,  0.000330643059,  0.000365280313, 
	 0.000403540267,  0.000445800921,  0.000492479594,  0.000544036797, 
	 0.000600981002,  0.000663873332,  0.000733333232,  0.000810044177, 
	  0.00089476048,   0.00098831451,   0.00109162449,   0.00120570394, 
	  0.00133167044,   0.00147075707,   0.00162432413,   0.00179387128, 
	   0.0019810528,   0.00218769279,   0.00241580233,   0.00266759819, 
	  0.00294552394,   0.00325227319,   0.00359081314,   0.00396441435, 
	  0.00437667919,   0.00483157486,   0.00533347158,   0.00588718057, 
	  0.00649799826,   0.00717175426,   0.00791486353,   0.00873438362, 
	  0.00963807665,      0.01063448,    0.0117329769,    0.0129438806, 
	   0.0142785236,    0.0157493539,    0.0173700359,    0.0191555731, 
	   0.0211224295,    0.0232886579,    0.0256740525,    0.0283003151, 
	   0.0311912112,    0.0343727656,    0.0378734544,     0.041724436, 
	    0.045959767,    0.0506166518,    0.0557357259,    0.0613613129, 
	   0.0675417557,    0.0743297189,    0.0817825496,    0.0899626091, 
	   0.0989376903,     0.108781397,     0.119573571,     0.131400704, 
	    0.144356415,     0.158541903,     0.174066424,     0.191047713, 
	    0.209612563,     0.229897171,     0.252047658,      0.27622056, 
	    0.302583039,     0.331313431,     0.362601429,     0.396648288, 
	    0.433666974,     0.473882198,     0.517530143,     0.564858317, 
	    0.616124868,     0.671597898,     0.731554449,     0.796279192, 
	     0.86606276,     0.941199422,      1.02198493,      1.10871315, 
	     1.20167255,      1.30114222,       1.4073869,      1.52065134, 
	     1.64115465,      1.76908231,       1.9045794,      2.04774189, 
	     2.19860673,      2.35714316,      2.52324057,      2.69669867, 
	     2.87721562,      3.06437612,      3.25763988,       3.4563303, 
	     3.65962362,      3.86654019,      4.07593489,      4.28649282, 
	     4.49672508,      4.70496798,      4.90938902,      5.10799122, 
	     5.29863167,      5.47903585,      5.64682627,      5.79955292, 
	     5.93473244,      6.04989433,      6.14263391,      6.21066809, 
	     6.25190353,      6.26450062,      6.24694729,      6.19812489, 
	     6.11738157,      6.00459146,      5.86021185,       5.6853261, 
	     5.48167372,      5.25166035,       4.9983511,      4.72543716, 
	     4.43717909,      4.13832712,      3.83401561,      3.52963567, 
	     3.23068929,      2.94262838,      2.67068624,      2.41970658, 
	     2.19397998,      1.99709332,      1.83180165,      1.69992697, 
	     1.60229325,      1.53869748,      1.50792205,      1.50778484, 
	     1.53522897,       1.5864414,       1.6570009,       1.7420435, 
	     1.83644021,      1.93497717,      2.03253412,      2.12425137, 
	     2.20568085,      2.27292037,      2.32271957,       2.3525641, 
	     2.36072969,       2.3463068,        2.309196,      2.25007486, 
	     2.17033482,      2.07199597,      1.95759714,      1.83007252, 
	     1.69261634,      1.54854429,      1.40115809,      1.25361896, 
	     1.10883725,     0.969381273,     0.837409914,     0.714629769, 
	    0.602277756,     0.501127362,      0.41151607,     0.333390117, 
	     0.26636216,     0.209777206,     0.162781686,     0.124391638, 
	   0.0935558602,    0.0692112222,    0.0503282957,    0.0359461084, 
	   0.0251962356,    0.0173167586,    0.0116575742,   0.00767859351, 
	  0.00494270399,    0.0031051999,   0.00190126954,   0.00113284099, 
	 0.000655786716,  0.000368195149,   0.00020013639,  0.000105117593, 
	 5.32419435e-05,  2.59506705e-05,  1.21452331e-05,  5.44541308e-06, 
	 2.33338528e-06,  9.53217352e-07,  3.70272346e-07,  1.36395798e-07, 
	 4.75125326e-08,  1.56051598e-08,  4.81784879e-09, 
/* Z= 7 -> N */
	              0,  6.60882943e-05,  7.30210159e-05,  8.06805256e-05, 
	 8.91429227e-05,  9.84922881e-05,  0.000108821492,   0.00012023309, 
	 0.000132840374,  0.000146768463,  0.000162155557,  0.000179154245, 
	 0.000197933085,  0.000218678237,  0.000241595204,  0.000266910996, 
	 0.000294876227,  0.000325767644,  0.000359890808,   0.00039758315, 
	  0.00043921711,  0.000485203898,  0.000535997504,  0.000592099095, 
	 0.000654061791,  0.000722496188,    0.0007980761,   0.00088154542, 
	  0.00097372476,   0.00107551995,   0.00118793023,    0.0013120583, 
	   0.0014491207,   0.00160045933,   0.00176755432,    0.0019520385, 
	  0.00215571234,   0.00238056132,   0.00262877462,   0.00290276599, 
	  0.00320519623,   0.00353899761,   0.00390740251,   0.00431397185, 
	  0.00476262951,   0.00525769824,   0.00580393802,   0.00640659127, 
	  0.00707143033,   0.00780480914,   0.00861372054,   0.00950585864, 
	   0.0104896883,    0.0115745207,    0.0127705932,    0.0140891569, 
	   0.0155425817,    0.0171444546,    0.0189096984,    0.0208547004, 
	   0.0229974426,    0.0253576618,    0.0279569961,    0.0308191776, 
	    0.033970207,     0.037438564,    0.0412554368,    0.0454549491, 
	   0.0500744246,    0.0551546589,    0.0607402287,     0.066879794, 
	   0.0736264437,    0.0810380653,    0.0891777128,    0.0981140435, 
	     0.10792172,      0.11868187,     0.130482584,       0.1434194, 
	    0.157595798,     0.173123762,     0.190124333,     0.208728135, 
	    0.229075983,     0.251319379,     0.275621086,     0.302155763, 
	    0.331110299,     0.362684369,     0.397090882,     0.434556156, 
	    0.475320399,     0.519637525,     0.567775488,     0.620015621, 
	    0.676652789,     0.737994134,      0.80435878,     0.876075923, 
	     0.95348376,       1.0369271,      1.12675488,      1.22331703, 
	     1.32696104,      1.43802726,      1.55684435,      1.68372262, 
	     1.81894827,      1.96277571,      2.11541867,       2.2770412, 
	     2.44774771,      2.62757134,      2.81646299,      3.01427722, 
	     3.22076058,      3.43553662,       3.6580925,      3.88776588, 
	     4.12373018,      4.36498308,      4.61033583,      4.85840225, 
	     5.10759497,      5.35611916,      5.60197544,      5.84296417, 
	      6.0766964,       6.3006115,      6.51200199,      6.70804548, 
	     6.88584232,      7.04246807,      7.17502737,       7.2807188, 
	     7.35690832,      7.40120554,      7.41154814,      7.38628578, 
	       7.324265,      7.22491074,      7.08830023,      6.91522884, 
	     6.70725679,      6.46674109,      6.19684267,      5.90150881, 
	     5.58542442,      5.25393772,      4.91295099,      4.56878138, 
	     4.22800016,      3.89724422,      3.58301234,      3.29145575, 
	     3.02816629,      2.79797482,       2.6047709,      2.45135117, 
	     2.33930421,      2.26894355,      2.23928571,      2.24808192, 
	      2.2918992,      2.36624861,      2.46575189,      2.58434439, 
	     2.71549582,      2.85244918,      2.98845792,      3.11701822, 
	     3.23208642,      3.32826996,      3.40099239,      3.44662142, 
	     3.46256042,      3.44729924,       3.4004221,      3.32257795, 
	     3.21540737,      3.08143735,      2.92394423,       2.7467916, 
	     2.55425501,      2.35083652,      2.14108086,      1.92940557, 
	       1.719944,      1.51641643,      1.32202828,      1.13940179, 
	    0.970538676,     0.816815376,     0.679006398,     0.557332039, 
	    0.451525003,     0.360910594,     0.284494013,     0.221049413, 
	    0.169205651,     0.127524376,    0.0945676938,    0.0689531937, 
	    0.049395632,    0.0347354971,     0.023955375,    0.0161858089, 
	   0.0107025569,   0.00691741239,    0.0043646032,   0.00268468633, 
	  0.00160750956,  0.000935512362,  0.000528284174,  0.000288972195, 
	 0.000152835812,  7.80094852e-05,  3.83492807e-05,  1.81196647e-05, 
	 8.21072899e-06,   3.5601347e-06,  1.47358628e-06,  5.80812127e-07, 
	 2.17433396e-07,  7.71041258e-08,  2.58266972e-08,  8.14752976e-09, 
	 2.41334908e-09,  6.69062583e-10,  1.73031353e-10, 
/* Z= 8 -> O */
	              0,  0.000100438039,  0.000110967849,  0.000122600773, 
	 0.000135452254,  0.000149649757,  0.000165334102,  0.000182660791, 
	 0.000201801566,  0.000222946066,  0.000246303709,  0.000272105826, 
	 0.000300607731,   0.00033209144,  0.000366868335,  0.000405282161, 
	  0.00044771246,  0.000494578329,  0.000546342286,  0.000603515131, 
	 0.000666660431,  0.000736400601,  0.000813422317,  0.000898483617, 
	 0.000992421061,   0.00109615782,   0.00121071259,   0.00133720948, 
	  0.00147688854,   0.00163111801,    0.0018014072,   0.00198942074, 
	  0.00219699508,   0.00242615491,   0.00267913355,   0.00295839296, 
	  0.00326664769,    0.0036068901,   0.00398241822,   0.00439686654, 
	  0.00485424139,   0.00535895443,   0.00591586763,   0.00653033471, 
	  0.00720825186,    0.0079561118,   0.00878105965,   0.00969096273, 
	   0.0106944749,     0.011801118,    0.0130213648,     0.014366732, 
	    0.015849879,    0.0174847208,    0.0192865431,    0.0212721396, 
	   0.0234599467,     0.025870204,    0.0285251196,    0.0314490497, 
	   0.0346687064,    0.0382133573,    0.0421150699,    0.0464089587, 
	   0.0511334501,    0.0563305803,    0.0620463081,    0.0683308393, 
	   0.0752390102,    0.0828306526,    0.0911710113,     0.100331172, 
	    0.110388562,      0.12142738,     0.133539185,     0.146823347, 
	    0.161387727,     0.177349165,     0.194834173,     0.213979483, 
	     0.23493278,      0.25785327,     0.282912374,       0.3102943, 
	    0.340196848,     0.372831792,      0.40842548,     0.447219461, 
	     0.48947069,     0.535452068,     0.585452378,     0.639776468, 
	    0.698745012,      0.76269418,     0.831974924,     0.906952024, 
	    0.988002658,      1.07551479,      1.16988504,      1.27151525, 
	     1.38081014,      1.49817228,      1.62399805,      1.75867152, 
	     1.90255833,      2.05599809,      2.21929622,      2.39271402, 
	     2.57645893,      2.77067304,      2.97541928,      3.19066978, 
	     3.41629004,      3.65202451,      3.89748049,      4.15211248, 
	     4.41520596,      4.68586254,      4.96298361,      5.24525833, 
	      5.5311532,      5.81890154,       6.1065011,      6.39171124, 
	     6.67206049,      6.94485474,       7.2071991,      7.45602179, 
	     7.68810987,      7.90015459,      8.08880424,      8.25073147, 
	     8.38270092,      8.48165703,       8.5448122,      8.56974125, 
	     8.55448341,       8.4976387,      8.39846706,      8.25697422, 
	     8.07399559,      7.85125065,      7.59138775,      7.29799652, 
	     6.97559023,      6.62955904,      6.26608181,      5.89200449, 
	     5.51468229,       5.1417861,      4.78108644,      4.44020844, 
	     4.12638187,      3.84618521,      3.60529971,      3.40828705, 
	     3.25839925,      3.15743375,      3.10564089,      3.10169101, 
	     3.14270592,      3.22435069,      3.34098577,      3.48587084, 
	     3.65141225,      3.82944059,      4.01150846,      4.18919277, 
	     4.35438967,       4.4995923,      4.61813831,      4.70441818, 
	     4.75404167,       4.7639513,      4.73248339,      4.65937567, 
	     4.54571962,      4.39386368,      4.20727205,      3.99034071, 
	     3.74819088,      3.48643827,      3.21095872,      2.92765641, 
	     2.64224625,      2.36006308,      2.08589935,      1.82388318, 
	     1.57739449,      1.34902275,      1.14056265,     0.953044772, 
	    0.786796212,     0.641522825,     0.516408086,      0.41022101, 
	    0.321425438,     0.248286113,     0.188966006,     0.141611308, 
	    0.104421958,    0.0757064298,    0.0539208874,    0.0376936644, 
	   0.0258366689,    0.0173458848,    0.0113932518,   0.00731229875, 
	  0.00457976153,   0.00279515958,   0.00165997457,  0.000957741751, 
	 0.000535960949,  0.000290406169,  0.000152083085,  7.68309692e-05, 
	 3.73690964e-05,  1.74627148e-05,  7.82339703e-06,  3.35261052e-06, 
	 1.37104132e-06,  5.33736568e-07,  1.97284677e-07,  6.90526107e-08, 
	 2.28225971e-08,  7.10185599e-09,  2.07428497e-09,  5.66847624e-10, 
	 1.44450993e-10,   3.4208171e-11,  7.50134416e-12, 
/* Z= 9 -> F */
	              0,  0.000145245256,  0.000160462718,  0.000177273207, 
	  0.00019584324,   0.00021635677,  0.000239016896,  0.000264047936, 
	 0.000291697506,  0.000322239153,    0.0003559748,   0.00039323789, 
	 0.000434396527,  0.000479857117,  0.000530068413,  0.000585525704, 
	 0.000646775821,  0.000714422378,  0.000789131562,  0.000871638767, 
	 0.000962755526,    0.0010633776,   0.00117449334,    0.0012971937, 
	  0.00143268181,   0.00158228562,   0.00174747012,   0.00192985113, 
	  0.00213121134,   0.00235351641,   0.00259893527,   0.00286985817, 
	  0.00316892215,    0.0034990334,   0.00386339659,   0.00426554354, 
	  0.00470936671,    0.0051991553,   0.00573963625,   0.00633601611, 
	  0.00699402951,   0.00771999313,   0.00852086022,   0.00940428954, 
	   0.0103787063,    0.0114533864,    0.0126385307,    0.0139453681, 
	   0.0153862387,    0.0169747155,    0.0187257156,     0.020655632, 
	   0.0227824729,    0.0251260139,    0.0277079679,    0.0305521619, 
	    0.033684738,     0.037134368,    0.0409324728,    0.0451134965, 
	   0.0497151501,    0.0547787249,    0.0603494048,    0.0664765909, 
	   0.0732142925,    0.0806214884,    0.0887625813,     0.097707808, 
	     0.10753376,     0.118323848,      0.13016887,      0.14316757, 
	    0.157427222,     0.173064291,     0.190205038,     0.208986282, 
	    0.229555964,     0.252074063,     0.276713073,     0.303658992, 
	    0.333111912,     0.365286767,      0.40041405,     0.438740402, 
	    0.480529338,     0.526061654,      0.57563585,     0.629568577, 
	    0.688194573,     0.751866698,     0.820955753,     0.895849705, 
	     0.97695303,      1.06468523,      1.15947914,      1.26177919, 
	     1.37203801,      1.49071336,      1.61826396,      1.75514472, 
	     1.90180039,      2.05865979,      2.22612715,      2.40457416, 
	     2.59432912,      2.79566669,      3.00879455,      3.23384118, 
	     3.47083998,       3.7197144,      3.98026109,      4.25213242, 
	      4.5348177,      4.82762718,      5.12967205,      5.43984747, 
	     5.75681829,      6.07900238,      6.40456152,      6.73139238, 
	     7.05712414,       7.3791194,      7.69448423,      8.00008392, 
	     8.29256725,       8.5684042,      8.82392502,      9.05538273, 
	     9.25901508,      9.43112564,      9.56817055,      9.66685963, 
	     9.72425938,      9.73790741,      9.70592499,      9.62713242, 
	     9.50115204,      9.32851124,      9.11071682,      8.85031986, 
	     8.55094242,      8.21727943,      7.85506248,      7.47097731, 
	     7.07254505,      6.66795826,      6.26587486,      5.87517738, 
	     5.50470257,      5.16295004,      4.85778189,      4.59612799, 
	     4.38370705,      4.22478294,      4.12196636,      4.07607317, 
	     4.08605337,      4.14899063,      4.26018095,      4.41328096, 
	     4.60052919,      4.81302309,       5.0410409,       5.2743969, 
	     5.50281096,      5.71627617,       5.9054122,      6.06178427, 
	      6.1781826,      6.24884701,      6.26963043,      6.23809624, 
	     6.15354872,      6.01699257,      5.83102989,      5.59969759, 
	      5.3282547,      5.02292967,      4.69064283,      4.33871651, 
	      3.9745841,      3.60551882,      3.23838711,      2.87943935, 
	     2.53414536,      2.20707822,      1.90184724,      1.62107921, 
	     1.36644125,       1.1387037,     0.937830389,      0.76309222, 
	    0.613192797,     0.486400068,     0.380675972,     0.293797851, 
	    0.223467767,     0.167405516,     0.123424374,    0.0894883648, 
	   0.0637520552,    0.0445840396,    0.0305763241,    0.0205419566, 
	   0.0135034733,   0.00867480226,   0.00543893501,   0.00332357036, 
	  0.00197649584,   0.00114213396,  0.000640272745,  0.000347618188, 
	 0.000182455769,  9.24108026e-05,  4.50768275e-05,  2.11332499e-05, 
	  9.5025016e-06,  4.08886081e-06,  1.67976589e-06,    6.572335e-07, 
	 2.44291556e-07,  8.60314273e-08,  2.86255855e-08,  8.97292374e-09, 
	 2.64163669e-09,  7.28101412e-10,   1.8726333e-10,  4.47879615e-11, 
	 9.92585181e-12,   2.0308102e-12,  3.82131935e-13, 
/* Z=10 -> Ne */
	              0,  0.000202054522,   0.00022320899,  0.000246576208, 
	 0.000272387289,   0.00030089746,  0.000332388532,  0.000367171597, 
	 0.000405590254,  0.000448023755,  0.000494890846,  0.000546653871, 
	 0.000603823166,  0.000666962063,  0.000736692513,  0.000813700841, 
	 0.000898744678,  0.000992660294,   0.00109637028,   0.00121089327, 
	  0.00133735279,   0.00147698878,   0.00163116946,   0.00180140405, 
	  0.00198935741,   0.00219686655,   0.00242595701,   0.00267886301, 
	   0.0029580486,   0.00326623046,   0.00360640418,     0.003981872, 
	  0.00439627422,   0.00485362299,   0.00535834022,   0.00591529766, 
	  0.00652986346,   0.00720795011,   0.00795607176,   0.00878139865, 
	  0.00969182793,    0.0106960526,    0.0118036401,     0.013025118, 
	   0.0143720694,    0.0158572327,     0.017494617,    0.0192996245, 
	   0.0212891828,    0.0234818906,    0.0258981809,    0.0285604857, 
	   0.0314934328,     0.034724053,    0.0382819921,    0.0421997607, 
	   0.0465129949,    0.0512607358,    0.0564857461,    0.0622348227, 
	   0.0685591772,    0.0755148008,    0.0831628889,    0.0915702656, 
	     0.10080988,     0.110961311,     0.122111283,     0.134354264, 
	    0.147793025,     0.162539363,     0.178714722,     0.196450844, 
	    0.215890646,     0.237188846,     0.260512829,     0.286043406, 
	    0.313975662,     0.344519794,     0.377901882,     0.414364755, 
	    0.454168707,     0.497592241,     0.544932723,     0.596506953, 
	    0.652651727,     0.713723898,     0.780100763,     0.852179885, 
	    0.930378556,      1.01513338,      1.10689914,      1.20614755, 
	     1.31336439,       1.4290483,      1.55370653,      1.68785167, 
	     1.83199644,       1.9866488,      2.15230417,      2.32943916, 
	     2.51850152,      2.71990132,      2.93399858,      3.16109109, 
	      3.4014008,      3.65505815,      3.92208552,      4.20238018, 
	     4.49569416,      4.80161667,      5.11955214,      5.44870186, 
	     5.78804255,      6.13630772,      6.49197102,      6.85323048, 
	     7.21799612,      7.58388329,      7.94820929,       8.3079977, 
	     8.65998745,      9.00065899,      9.32625771,       9.6328392, 
	     9.91631794,      10.1725388,      10.3973494,      10.5866909, 
	     10.7367058,      10.8438435,      10.9049826,      10.9175673, 
	     10.8797178,      10.7903757,      10.6494122,      10.4577332, 
	     10.2173738,      9.93154907,      9.60468674,      9.24241734, 
	     8.85151577,      8.43980312,      8.01599693,      7.58950949, 
	     7.17020178,      6.76809978,      6.39307594,      6.05450821, 
	     5.76093626,      5.51972151,      5.33673525,      5.21608829, 
	     5.15991306,      5.16822386,      5.23885155,      5.36746931, 
	     5.54770803,      5.77135468,      6.02863455,      6.30855513, 
	     6.59930372,      6.88867712,      7.16452456,      7.41518593, 
	     7.62990475,       7.7991991,      7.91517496,      7.97177505, 
	     7.96494532,      7.89272499,      7.75524473,      7.55465269, 
	     7.29495335,       6.9817872,       6.6221509,      6.22407532, 
	     5.79628325,      5.34783602,      4.88779116,      4.42488527, 
	     3.96725726,      3.52221465,      3.09606194,      2.69398451, 
	     2.31999016,      1.97690785,      1.66643286,       1.3892138, 
	     1.14497101,     0.932635844,     0.750502765,     0.596383929, 
	     0.46775949,     0.361917019,     0.276075065,     0.207487971, 
	    0.153530374,     0.111761093,    0.0799670815,    0.0561894588, 
	   0.0387337096,    0.0261667743,    0.0173038561,    0.0111876931, 
	  0.00706293201,    0.0043479749,   0.00260630855,   0.00151896977, 
	 0.000859355379,  0.000471172563,  0.000249931123,  0.000128029104, 
	 6.32155279e-05,  3.00266001e-05,  1.36916751e-05,  5.98044971e-06, 
	 2.49662025e-06,  9.93757794e-07,  3.76216178e-07,  1.35111264e-07, 
	 4.59048337e-08,  1.47129171e-08,  4.43524062e-09,  1.25360589e-09, 
	 3.31146305e-10,  8.14753195e-11,  1.86061114e-11,   3.9294349e-12, 
	 7.64573843e-13,  1.36536805e-13,  2.22895382e-14, 
/* Z=11 -> Na */
	              0,  0.000273550104,  0.000302168075,  0.000333776901, 
	 0.000368688692,  0.000407247979,  0.000449835206,  0.000496870372, 
	 0.000548817101,   0.00060618727,  0.000669545901,  0.000739516632, 
	 0.000816787884,    0.0009021195,  0.000996350078,   0.00110040512, 
	  0.00121530588,    0.0013421796,   0.00148226961,    0.0016369482, 
	  0.00180772913,   0.00199628249,   0.00220445031,    0.0024342651, 
	  0.00268796762,   0.00296802958,   0.00327717676,   0.00361841451, 
	   0.0039950558,   0.00441075349,    0.0048695337,   0.00537583418, 
	  0.00593454437,   0.00655105198,   0.00723129371,   0.00798180699, 
	  0.00880979374,    0.0097231837,     0.010730708,    0.0118419761, 
	   0.0130675649,    0.0144191142,    0.0159094259,    0.0175525844, 
	   0.0193640757,    0.0213609189,    0.0235618278,    0.0259873588, 
	   0.0286600925,    0.0316048227,    0.0348487683,    0.0384217948, 
	   0.0423566736,    0.0466893241,    0.0514591262,    0.0567092225, 
	   0.0624868646,     0.068843767,    0.0758365095,    0.0835269615, 
	   0.0919827372,     0.101277687,     0.111492403,      0.12271481, 
	    0.135040745,     0.148574561,     0.163429856,     0.179730132, 
	    0.197609574,     0.217213795,     0.238700673,     0.262241215, 
	    0.288020372,     0.316237986,     0.347109675,      0.38086772, 
	    0.417762071,     0.458061099,     0.502052605,     0.550044417, 
	    0.602365434,     0.659365952,     0.721418381,     0.788917542, 
	    0.862280846,     0.941948056,      1.02838111,       1.1220634, 
	     1.22349834,      1.33320808,      1.45173144,       1.5796206, 
	     1.71743822,      1.86575282,      2.02513361,      2.19614458, 
	     2.37933707,      2.57524109,      2.78435588,      3.00713944, 
	     3.24399471,      3.49525762,      3.76117921,       4.0419116, 
	     4.33748722,       4.6478014,      4.97258997,        5.311409, 
	     5.66361094,      6.02832508,      6.40443182,      6.79054499, 
	     7.18499088,      7.58579159,      7.99065065,      8.39694691, 
	      8.8017292,      9.20172215,      9.59333706,       9.9726963, 
	     10.3356667,      10.6779022,      10.9949064,      11.2821007, 
	     11.5349188,      11.7488966,      11.9197979,      12.0437307, 
	     12.1172915,      12.1376963,      12.1029377,      12.0119171, 
	      11.864583,      11.6620464,      11.4066811,      11.1021872, 
	     10.7536249,      10.3673983,      9.95119381,      9.51386547, 
	     9.06525898,      8.61598587,      8.17713833,      7.75996208, 
	     7.37548876,      7.03414583,      6.74535513,      6.51714325, 
	     6.35577726,      6.26544905,      6.24802685,      6.30288935, 
	     6.42685223,       6.6142025,      6.85683632,      7.14450169, 
	     7.46513939,      7.80529976,      8.15063095,      8.48640156, 
	     8.79804993,       9.0717268,      9.29480076,      9.45632458, 
	     9.54742527,      9.56160831,      9.49495888,      9.34624863, 
	      9.1169157,      8.81095123,      8.43467903,      7.99644518, 
	     7.50623369,      6.97523212,       6.4153614,        5.838799, 
	       5.257514,      4.68284035,      4.12510538,       3.5933249, 
	     3.09498072,      2.63588238,      2.22011399,      1.85006249, 
	     1.52651763,      1.24883115,      1.01512277,     0.822514594, 
	    0.667379439,     0.545587838,     0.452740431,     0.384374678, 
	    0.336139709,      0.30393374,     0.284004658,     0.273014754, 
	    0.268074125,     0.266748607,     0.267047912,     0.267400444, 
	    0.266619921,     0.263868004,     0.258615553,     0.250604331, 
	    0.239809483,     0.226402283,     0.210712954,     0.193192944, 
	    0.174376294,     0.154841036,     0.135171324,     0.115922071, 
	   0.0975879431,    0.0805787295,     0.065202266,    0.0516563915, 
	    0.040029645,    0.0303101279,    0.0224010423,    0.0161406845, 
	   0.0113245118,   0.00772683416,   0.00512009906,   0.00329023856, 
	  0.00204735063,   0.00123163359,  0.000715099042,  0.000400017307, 
	 0.000215184511,  0.000111099413,  5.49395663e-05, 
/* Z=12 -> Mg */
	              0,  0.000361421058,  0.000399200799,  0.000440925272, 
	 0.000487005687,  0.000537895947,   0.00059409719,  0.000656162505, 
	 0.000724702433,  0.000800391019,  0.000883971981,  0.000976266398, 
	  0.00107818015,   0.00119071337,   0.00131496915,   0.00145216507, 
	  0.00160364457,   0.00177088974,    0.0019555355,   0.00215938571, 
	  0.00238443026,   0.00263286335,   0.00290710549,   0.00320982584, 
	  0.00354396761,   0.00391277578,   0.00431982847,   0.00476906914, 
	  0.00526484521,   0.00581194693,   0.00641565397,   0.00708178338, 
	  0.00781674311,   0.00862759165,    0.0095221037,     0.010508839, 
	   0.0115972208,    0.0127976257,      0.01412147,    0.0155813163, 
	   0.0171909854,    0.0189656802,    0.0209221132,    0.0230786614, 
	   0.0254555214,    0.0280748829,    0.0309611186,    0.0341409966, 
	   0.0376438983,     0.041502066,    0.0457508676,    0.0504290946, 
	   0.0555792563,    0.0612479337,    0.0674861446,    0.0743497312, 
	   0.0818997994,    0.0902031511,    0.0993328243,      0.10936857, 
	    0.120397471,     0.132514492,     0.145823181,     0.160436332, 
	    0.176476687,     0.194077775,     0.213384613,     0.234554693, 
	    0.257758737,     0.283181727,     0.311023802,     0.341501206, 
	    0.374847353,     0.411313802,     0.451171219,     0.494710475, 
	    0.542243481,     0.594104111,     0.650649071,     0.712258637, 
	    0.779337287,     0.852313936,     0.931642413,      1.01780128, 
	     1.11129355,      1.21264589,       1.3224076,      1.44114888, 
	     1.56945872,      1.70794189,      1.85721529,      2.01790357, 
	     2.19063354,      2.37602758,       2.5746963,      2.78722906, 
	       3.014184,       3.2560761,      3.51336408,      3.78643489, 
	     4.07558775,      4.38101625,      4.70278788,      5.04082251, 
	     5.39487171,      5.76449299,      6.14902782,      6.54757452, 
	     6.95896673,      7.38174677,      7.81414795,      8.25407124, 
	     8.69907379,      9.14635372,      9.59275055,      10.0347452, 
	     10.4684715,      10.8897448,      11.2940874,      11.6767817, 
	     12.0329342,      12.3575478,      12.6456194,      12.8922462, 
	     13.0927525,      13.2428217,      13.3386536,      13.3771143, 
	      13.355896,      13.2736826,      13.1302938,      12.9268217, 
	     12.6657476,       12.351017,      11.9880829,      11.5839014, 
	     11.1468668,      10.6866951,      10.2142305,      9.74120808, 
	     9.27993107,      8.84291267,      8.44246101,      8.09023762, 
	     7.79679918,      7.57114458,      7.42029238,      7.34890556, 
	     7.35898876,      7.44968081,      7.61715174,      7.85462236, 
	     8.15251064,      8.49870014,      8.87892532,      9.27726746, 
	     9.67672157,      10.0598316,      10.4093552,      10.7089233, 
	     10.9436798,      11.1008606,        11.17029,      11.1447639, 
	     11.0203276,      10.7963982,      10.4757557,       10.064395, 
	     9.57123184,      9.00769901,      8.38723564,      7.72470713, 
	     7.03576994,      6.33623028,      5.64141607,      4.96559954, 
	     4.32149315,      3.71984911,      3.16916871,      2.67553949, 
	       2.242594,      1.87159145,      1.56160271,       1.3097831, 
	     1.11170924,     0.961754739,     0.853479385,     0.780007541, 
	     0.73437655,     0.709839046,     0.700109482,     0.699550748, 
	    0.703300476,     0.707342148,     0.708527386,     0.704558551, 
	    0.693938971,     0.675899565,     0.650309324,     0.617573857, 
	    0.578528106,     0.534326911,     0.486336768,     0.436031997, 
	     0.38489908,     0.334352016,     0.285661131,     0.239898384, 
	    0.197900355,     0.160249591,     0.127273902,    0.0990620181, 
	   0.0754928365,    0.0562745184,    0.0409895629,    0.0291412994, 
	   0.0201979559,    0.0136309844,   0.00894531142,   0.00570042385, 
	  0.00352228107,   0.00210704259,   0.00121827784,   0.00067966833, 
	 0.000365207874,  0.000188646765,  9.34879281e-05,  4.43552162e-05, 
	 2.01027633e-05,  8.68317966e-06,  3.56576311e-06, 
/* Z=13 -> Al */
	              0,  0.000467162929,  0.000515953114,  0.000569832744, 
	 0.000629331684,   0.00069503492,  0.000767588033,  0.000847703544, 
	 0.000936167955,    0.0010338492,   0.00114170508,   0.00126079272, 
	  0.00139227835,   0.00153744884,   0.00169772422,   0.00187467085, 
	  0.00207001693,   0.00228566886,   0.00252372934,   0.00278651761, 
	  0.00307659106,   0.00339677045,   0.00375016569,   0.00414020615, 
	  0.00457067136,   0.00504572922,    0.0055699721,   0.00614846172, 
	  0.00678677624,     0.007491061,   0.00826808624,   0.00912530906, 
	   0.0100709433,     0.011114032,     0.012264532,    0.0135334032, 
	   0.0149327088,    0.0164757203,    0.0181770436,    0.0200527329, 
	   0.0221204534,    0.0243996195,    0.0269115698,    0.0296797529, 
	   0.0327299275,     0.036090374,    0.0397921503,    0.0438693315, 
	   0.0483593009,    0.0533030592,    0.0587455519,     0.064736031, 
	   0.0713284314,      0.07858181,    0.0865607932,    0.0953360572, 
	     0.10498485,     0.115591571,     0.127248377,     0.140055791, 
	    0.154123455,     0.169570774,     0.186527804,     0.205135971, 
	    0.225549012,      0.24793385,     0.272471547,     0.299358338, 
	    0.328806609,     0.361045927,     0.396324247,     0.434908867, 
	    0.477087557,     0.523169756,     0.573487401,     0.628396213, 
	     0.68827635,     0.753533483,     0.824599385,     0.901932478, 
	    0.986018181,      1.07736886,      1.17652404,      1.28404915, 
	     1.40053475,      1.52659512,      1.66286576,      1.81000078, 
	     1.96866918,      2.13954997,      2.32332706,      2.52068233, 
	     2.73228765,      2.95879579,      3.20082974,      3.45897031, 
	       3.733742,      4.02559757,      4.33490133,      4.66190863, 
	     5.00674582,       5.3693881,      5.74963379,      6.14708233, 
	     6.56110382,      6.99081564,      7.43505478,      7.89235163, 
	      8.3609066,      8.83856487,       9.3228054,      9.81071758, 
	     10.2990036,      10.7839708,      11.2615499,      11.7273064, 
	     12.1764812,      12.6040363,       13.004715,       13.373126, 
	     13.7038336,      13.9914789,      14.2309046,       14.417305, 
	     14.5463829,       14.614522,      14.6189623,      14.5579739, 
	     14.4310284,      14.2389526,      13.9840641,       13.670269, 
	     13.3031235,       12.889842,      12.4392509,      11.9616699, 
	     11.4687309,      10.9731178,      10.4882441,      10.0278511, 
	     9.60557175,      9.23442936,      8.92633152,      8.69155025, 
	     8.53823376,       8.4719553,      8.49534893,      8.60783291, 
	     8.80545902,      9.08089066,      9.42353153,      9.81980038, 
	     10.2535419,      10.7065754,       11.159338,      11.5916281, 
	     11.9833794,      12.3154602,      12.5704517,      12.7333622, 
	     12.7922564,      12.7387571,      12.5684052,      12.2808514, 
	     11.8798637,      11.3731689,      10.7721052,      10.0911207, 
	     9.34714127,      8.55881691,      7.74572134,      6.92750549, 
	     6.12308359,      5.34987116,      4.62312603,      3.95542169, 
	     3.35627842,      2.83196664,       2.3854847,      2.01670551, 
	     1.72267115,      1.49800766,      1.33542681,      1.22627699, 
	     1.16110373,       1.1301862,      1.12402046,      1.13372588, 
	     1.15136445,      1.17016315,      1.18464684,      1.19068456, 
	     1.18546295,      1.16740179,      1.13602138,      1.09178042, 
	     1.03589344,     0.970140815,      0.89668113,     0.817872941, 
	    0.736114979,     0.653708458,     0.572746456,     0.495031953, 
	    0.422025144,     0.354818791,      0.29413861,     0.240365028, 
	    0.193571374,     0.153573051,     0.119982883,    0.0922675431, 
	   0.0698012263,    0.0519134626,    0.0379291624,    0.0271998998, 
	   0.0191264693,    0.0131733082,   0.00887588784,   0.00584246451, 
	  0.00375151099,   0.00234612287,   0.00142656185,  0.000841868867, 
	 0.000481269788,  0.000265988259,  0.000141830067,  7.28069135e-05, 
	  3.5901001e-05,  1.69654377e-05,  7.66488029e-06, 
/* Z=14 -> Si */
	              0,  0.000592864642,  0.000654725125,  0.000723031757, 
	 0.000798454974,  0.000881734595,  0.000973687158,   0.00107521378, 
	  0.00118730881,   0.00131106959,   0.00144770672,   0.00159855606, 
	  0.00176509155,   0.00194893905,   0.00215189205,   0.00237592915, 
	  0.00262323255,   0.00289620901,   0.00319751259,   0.00353007065, 
	     0.00389711,   0.00430218875,   0.00474923011,    0.0052425568, 
	  0.00578693533,   0.00638761744,   0.00705039082,   0.00778163131, 
	  0.00858836435,   0.00947832596,    0.0104600396,    0.0115428884, 
	   0.0127372034,    0.0140543552,    0.0155068627,    0.0171084981, 
	   0.0188744143,    0.0208212752,    0.0229674131,    0.0253329743, 
	   0.0279401094,    0.0308131557,    0.0339788496,    0.0374665558, 
	   0.0413085073,    0.0455400981,    0.0502001382,     0.055331219, 
	   0.0609800145,    0.0671976879,    0.0740402788,     0.081569165, 
	   0.0898514986,    0.0989607573,     0.108977258,     0.119988777, 
	    0.132091165,     0.145389035,     0.159996465,     0.176037803, 
	    0.193648472,     0.212975815,     0.234180033,     0.257435173, 
	    0.282930076,     0.310869545,     0.341475219,     0.374987036, 
	    0.411664069,     0.451785862,     0.495653629,      0.54359138, 
	    0.595947146,     0.653093874,     0.715430915,     0.783384502, 
	       0.857409,     0.937987328,      1.02563143,      1.12088263, 
	     1.22431135,      1.33651721,      1.45812714,      1.58979475, 
	     1.73219788,      1.88603616,      2.05202746,      2.23090267, 
	     2.42340136,      2.63026452,       2.8522265,      3.09000611, 
	      3.3442955,      3.61574793,      3.90496325,      4.21247196, 
	     4.53871775,      4.88403749,      5.24863863,      5.63257742, 
	     6.03573179,      6.45777559,      6.89815187,      7.35604143, 
	     7.83033562,      8.31960964,      8.82209015,      9.33563614, 
	      9.8577137,      10.3853798,      10.9152708,      11.4436007, 
	     11.9661646,      12.4783592,      12.9752064,      13.4514046, 
	     13.9013891,      14.3194027,      14.6995993,      15.0361576, 
	     15.3234148,      15.5560188,      15.7290955,      15.8384304, 
	     15.8806562,      15.8534498,      15.7557154,      15.5877638, 
	     15.3514719,      15.0504084,      14.6899137,      14.2771473, 
	     13.8210478,      13.3322506,      12.8229084,      12.3064432, 
	     11.7972202,      11.3101368,      10.8601522,      10.4617586, 
	      10.128418,      9.87198162,      9.70212555,      9.62582779, 
	     9.64691162,      9.76568794,      9.97872734,      10.2787743, 
	     10.6548223,      11.0923624,      11.5737944,      12.0790005, 
	      12.586051,      13.0720282,       13.513916,      13.8895435, 
	     14.1785059,      14.3630438,      14.4288311,      14.3656273, 
	     14.1677675,      13.8344498,      13.3698158,      12.7828035, 
	     12.0867825,       11.298975,       10.439703,      9.53148174, 
	     8.59801292,      7.66313267,       6.7497673,      5.87895155, 
	     5.06897068,      4.33466911,      3.68696117,      3.13257432, 
	     2.67402267,      2.30981112,      2.03484011,      1.84097946, 
	     1.71775663,      1.65311468,      1.63417959,      1.64799058, 
	     1.68215203,      1.72537243,      1.76787519,      1.80167127, 
	     1.82069683,      1.82082844,      1.79979384,      1.75699639, 
	      1.6932795,       1.6106497,      1.51198208,      1.40072465, 
	     1.28061903,      1.15544927,      1.02882969,     0.904036164, 
	    0.783884823,     0.670657516,     0.566069841,     0.471277267, 
	    0.386911988,     0.313142091,     0.249745756,     0.196192533, 
	    0.151725456,     0.115438528,    0.0863461569,    0.0634417832, 
	   0.0457450412,    0.0323370472,    0.0223847087,    0.0151551683, 
	   0.0100218868,   0.00646395609,   0.00406021671,   0.00247974182, 
	  0.00147007324,  0.000844462833,  0.000469170569,  0.000251625519, 
	 0.000130011729,  6.45815089e-05,  3.07742739e-05,  1.40358216e-05, 
	 6.11269206e-06,   2.5357208e-06,   9.9937165e-07, 
/* Z=15 -> P */
	              0,  0.000740636024,  0.000817838358,  0.000903076841, 
	 0.000997186173,   0.00110108731,    0.0012157968,   0.00134243595, 
	  0.00148224225,   0.00163658115,   0.00180695835,   0.00199503545, 
	  0.00220264564,   0.00243181037,   0.00268475921,   0.00296395156, 
	  0.00327209942,   0.00361219281,   0.00398752932,   0.00440174388, 
	  0.00485884398,   0.00536324689,   0.00591982156,   0.00653393334, 
	  0.00721149426,   0.00795901846,   0.00878368318,   0.00969339162, 
	   0.0106968507,    0.0118036484,    0.0130243357,    0.0143705346, 
	   0.0158550348,    0.0174919032,    0.0192966275,     0.021286238, 
	   0.0234794635,    0.0258968957,    0.0285611749,     0.031497173, 
	   0.0347322263,    0.0382963493,    0.0422225073,    0.0465468802, 
	   0.0513091721,    0.0565529391,    0.0623259284,     0.068680495, 
	     0.07567399,    0.0833692253,    0.0918349549,     0.101146422, 
	    0.111385897,     0.122643299,     0.135016873,      0.14861384, 
	    0.163551196,     0.179956526,     0.197968751,      0.21773918, 
	    0.239432365,     0.263227105,     0.289317608,     0.317914486, 
	    0.349245936,     0.383558989,     0.421120703,     0.462219387, 
	    0.507165909,     0.556294978,     0.609966457,     0.668566465, 
	     0.73250854,     0.802235067,     0.878217757,     0.960958898, 
	     1.05099177,      1.14888132,      1.25522363,      1.37064636, 
	     1.49580789,      1.63139582,      1.77812529,      1.93673682, 
	     2.10799241,      2.29267216,      2.49156785,      2.70547748, 
	     2.93519759,      3.18151331,      3.44518828,      3.72695303, 
	     4.02748871,      4.34741306,      4.68726158,      5.04746723, 
	     5.42833805,      5.83003378,      6.25253916,      6.69563723, 
	     7.15887737,      7.64154816,      8.14264393,      8.66083527, 
	     9.19443893,      9.74138832,      10.2992077,       10.864994, 
	     11.4353971,      12.0066147,      12.5743923,      13.1340313, 
	     13.6804171,      14.2080584,      14.7111387,      15.1835966, 
	     15.6192112,      16.0117188,      16.3549519,      16.6429825, 
	     16.8703041,      17.0320168,      17.1240215,      17.1432362, 
	     17.0878029,      16.9572792,      16.7528324,      16.4773846, 
	     16.1357441,       15.734664,      15.2828541,      14.7909174, 
	     14.2712078,       13.737606,      13.2051983,      12.6898746, 
	     12.2078514,      11.7751112,       11.406805,      11.1166153, 
	     10.9161148,      10.8141661,      10.8163614,      10.9245672, 
	     11.1365957,       11.446022,      11.8421707,      12.3103065, 
	     12.8320036,      13.3857069,      13.9474688,      14.4918318, 
	     14.9928207,      15.4250059,      15.7645874,      15.9904528, 
	      16.085146,      16.0357075,      15.8343391,      15.4788456, 
	     14.9728298,      14.3256168,       13.551918,       12.671217, 
	     11.7069273,      10.6853371,      9.63441849,       8.5825386, 
	     7.55715942,      6.58359528,      5.68390274,      4.87596273, 
	     4.17281675,      3.58228779,      3.10690594,       2.7441349, 
	     2.48686981,      2.32416701,      2.24213624,      2.22493482, 
	     2.25578451,       2.3179431,      2.39557457,      2.47446847, 
	     2.54258108,      2.59038544,      2.61103368,      2.60034466, 
	     2.55664372,      2.48048592,      2.37429166,      2.24193478, 
	     2.08830929,      1.91890645,      1.73942161,      1.55541492, 
	     1.37203085,      1.19378972,      1.02444625,     0.866917729, 
	    0.723270595,     0.594757378,     0.481892198,     0.384553939, 
	     0.30210492,     0.233515173,     0.177483857,     0.132551104, 
	   0.0971955359,    0.0699145645,    0.0492866151,    0.0340152122, 
	   0.0229564272,    0.0151317073,   0.00972862635,   0.00609235093, 
	  0.00371057494,   0.00219450681,   0.00125820772,  0.000698123302, 
	 0.000374183088,  0.000193365966,  9.61509068e-05,  4.59090916e-05, 
	 2.10024446e-05,  9.18497608e-06,  3.83080305e-06,  1.51993436e-06, 
	 5.72207568e-07,  2.03842717e-07,  6.85196682e-08, 
/* Z=16 -> S */
	              0,   0.00091274752,   0.00100779033,   0.00111271499, 
	  0.00122854637,    0.0013564158,   0.00149757077,   0.00165338814, 
	  0.00182538654,   0.00201524119,   0.00222480041,   0.00245610229, 
	  0.00271139457,    0.0029931569,   0.00330412271,    0.0036473074, 
	  0.00402603392,   0.00444396771,   0.00490514748,   0.00541402632, 
	  0.00597551186,   0.00659501273,   0.00727848848,    0.0080325054, 
	  0.00886429753,   0.00978183281,    0.0107938861,    0.0119101228, 
	   0.0131411841,    0.0144987796,    0.0159958061,    0.0176464487, 
	   0.0194663201,    0.0214725919,    0.0236841496,    0.0261217616, 
	   0.0288082585,    0.0317687318,    0.0350307524,    0.0386246033, 
	   0.0425835475,    0.0469440892,    0.0517463088,    0.0570341609, 
	   0.0628558695,    0.0692642853,    0.0763173401,    0.0840784833, 
	   0.0926171988,     0.102009527,      0.11233864,     0.123695493, 
	    0.136179447,     0.149899021,     0.164972648,     0.181529537, 
	    0.199710444,     0.219668716,     0.241571262,     0.265599489, 
	    0.291950524,     0.320838362,     0.352494895,     0.387171477, 
	    0.425139904,     0.466693878,     0.512150407,     0.561851144, 
	    0.616163611,      0.67548269,     0.740231991,     0.810864806, 
	    0.887865603,     0.971750736,      1.06306946,      1.16240442, 
	     1.27037215,      1.38762248,      1.51483941,      1.65273857, 
	     1.80206692,      1.96360004,      2.13813949,      2.32650852, 
	     2.52954793,      2.74810934,       2.9830482,      3.23521495, 
	       3.505445,       3.7945466,      4.10328627,      4.43237543, 
	     4.78244781,       5.1540451,      5.54759026,      5.96336555, 
	     6.40148449,      6.86186409,      7.34419537,       7.8479104, 
	     8.37214851,      8.91572666,      9.47710133,      10.0543451, 
	     10.6451054,      11.2465906,      11.8555431,      12.4682245, 
	     13.0804119,      13.6874027,        14.28403,      14.8646946, 
	     15.4234161,      15.9538908,      16.4495888,      16.9038544, 
	     17.3100357,      17.6616402,      17.9525051,      18.1769886, 
	     18.3301792,      18.4081135,      18.4079971,      18.3284321, 
	     18.1696262,      17.9335728,      17.6242199,      17.2475643, 
	     16.8117142,      16.3268623,      15.8051939,      15.2606869, 
	     14.7088366,      14.1662731,      13.6502857,      13.1782684, 
	     12.7670841,      12.4323969,      12.1879616,       12.044939, 
	     12.0112352,      12.0909376,      12.2838593,       12.585228, 
	     12.9855642,      13.4707537,      14.0223455,      14.6180544, 
	     15.2324743,      15.8379927,      16.4058475,      16.9072914, 
	     17.3148327,      17.6034641,      17.7518444,      17.7433586, 
	     17.5669937,      17.2179985,      16.6982479,      16.0163269, 
	     15.1872616,      14.2319593,      13.1763077,      12.0500317, 
	     10.8853197,      9.71530819,      8.57251072,      7.48727083, 
	     6.48633814,      5.59166098,      4.81946135,      4.17965841, 
	     3.67566919,      3.30459619,      3.05777526,      2.92163873, 
	     2.87881637,      2.90938902,       2.9921999,      3.10612512, 
	      3.2312274,      3.34971929,      3.44669271,      3.51059222, 
	     3.53343177,      3.51076841,      3.44147086,      3.32732391, 
	     3.17251444,      2.98305154,      2.76616383,      2.52971792, 
	      2.2816875,      2.02970147,      1.78068399,      1.54059505, 
	     1.31426895,       1.1053462,     0.916284442,     0.748433888, 
	    0.602161884,     0.477008641,     0.371858925,      0.28511709, 
	    0.214873254,     0.159053236,     0.115546107,    0.0823070109, 
	   0.0574336685,     0.039218165,     0.026176095,    0.0170564465, 
	   0.0108360397,   0.00670265267,   0.00403066538,    0.0023527944, 
	   0.0013309276,  0.000728349201,  0.000384903193,  0.000196049819, 
	 9.60551552e-05,  4.51762317e-05,  2.03511081e-05,  8.76121703e-06, 
	 3.59585897e-06,  1.40351801e-06,  5.19607454e-07,  1.81964765e-07, 
	 6.01051937e-08,  1.86703115e-08,  5.43689627e-09, 
/* Z=17 -> Cl */
	              0,   0.00111159275,   0.00122721237,   0.00135483884, 
	  0.00149571558,   0.00165121513,   0.00182285078,   0.00201229239, 
	  0.00222138152,   0.00245214929,   0.00270683598,   0.00298791216, 
	  0.00329810171,   0.00364040839,   0.00401814422,   0.00443496043, 
	  0.00489488151,   0.00540234428,   0.00596223772,   0.00657995045, 
	   0.0072614206,   0.00801318977,   0.00884246547,   0.00975718629, 
	   0.0107660964,    0.0118788239,    0.0131059708,    0.0144592104, 
	   0.0159513876,    0.0175966416,    0.0194105301,    0.0214101691, 
	   0.0236143805,    0.0260438677,    0.0287213884,    0.0316719599, 
	   0.0349230692,    0.0385049246,    0.0424506888,    0.0467967838, 
	   0.0515831858,    0.0568537675,    0.0626566485,    0.0690445825, 
	   0.0760754123,     0.083812505,    0.0923252702,     0.101689681, 
	     0.11198888,      0.12331377,     0.135763749,     0.149447396, 
	    0.164483249,     0.181000665,     0.199140653,     0.219056934, 
	    0.240916833,     0.264902383,     0.291211486,     0.320059091, 
	    0.351678342,     0.386322051,     0.424263865,     0.465799779, 
	    0.511249483,     0.560957849,     0.615296304,     0.674664497, 
	    0.739491463,      0.81023705,     0.887393355,     0.971485853, 
	     1.06307423,      1.16275334,      1.27115428,      1.38894355, 
	     1.51682365,      1.65553296,      1.80584347,      1.96855986, 
	     2.14451694,      2.33457589,      2.53962016,      2.76054978, 
	     2.99827528,       3.2537086,      3.52775407,      3.82129669, 
	     4.13518953,      4.47023678,        4.827178,      5.20666838, 
	     5.60925627,      6.03535843,      6.48523378,      6.95895672, 
	     7.45638227,      7.97711706,      8.52048206,      9.08548164, 
	     9.67076302,      10.2745867,      10.8947906,        11.52876, 
	     12.1733999,       12.825119,      13.4798126,       14.132863, 
	     14.7791424,      15.4130392,      16.0284977,      16.6190681, 
	     17.1779919,      17.6982918,      18.1728954,      18.5947876, 
	     18.9571648,      19.2536373,      19.4784374,      19.6266384, 
	     19.6944008,       19.679203,      19.5800838,      19.3978558, 
	     19.1352997,      18.7973251,      18.3910561,      17.9258785, 
	     17.4133797,      16.8672295,      16.3029346,      15.7375021, 
	     15.1890001,      14.6760092,      14.2169981,      13.8296156, 
	     13.5299473,      13.3317442,      13.2456808,      13.2786608, 
	     13.4332314,      13.7071228,      14.0929823,      14.5782957, 
	     15.1455545,      15.7726536,      16.4335423,      17.0990982, 
	     17.7382011,      18.3189812,      18.8101635,      19.1824818, 
	     19.4100513,      19.4716816,      19.3520069,      19.0424023, 
	     18.5416145,      17.8560371,      16.9996529,      15.9935503, 
	     14.8650942,      13.6467266,       12.374485,       11.086277, 
	     9.82003975,      8.61186886,       7.4942317,      6.49439001, 
	     5.63312721,      4.92386723,      4.37224722,      3.97616935, 
	     3.72631812,      3.60710049,      3.59792686,      3.67472887, 
	     3.81159472,      3.98239851,      4.16230583,      4.32906342, 
	     4.46399832,      4.55268478,      4.58527422,      4.55649281, 
	     4.46536398,       4.3146925,      4.11038971,      3.86070514, 
	     3.57542229,      3.26508594,      2.94029999,      2.61113477, 
	     2.28666115,      1.97462106,      1.68123245,      1.41111565, 
	     1.16732478,     0.951462448,       0.7638551,     0.603766382, 
	    0.469628483,     0.359273106,     0.270149827,     0.199519709, 
	    0.144619003,     0.102788813,     0.071570605,    0.0487689488, 
	    0.032485377,    0.0211275276,    0.0133989165,   0.00827483274, 
	  0.00496924715,   0.00289734756,   0.00163752632,  0.000895611709, 
	 0.000473172491,  0.000241033966,  0.000118153213,  5.56198902e-05, 
	 2.50898647e-05,  1.08211161e-05,  4.45171054e-06,    1.742569e-06, 
	 6.47348088e-07,  2.27611622e-07,  7.55316591e-08,  2.35860327e-08, 
	 6.90920077e-09,   1.8925419e-09,  4.83111384e-10, 
/* Z=18 -> Ar */
	              0,   0.00133971556,   0.00147889985,   0.00163251965, 
	   0.0018020682,   0.00198919256,   0.00219570962,   0.00242362265, 
	  0.00267514121,   0.00295270234,   0.00325899269,   0.00359697454, 
	  0.00396991428,   0.00438141311,   0.00483543938,   0.00533636892, 
	  0.00588902319,   0.00649871631,    0.0071713035,   0.00791323651, 
	  0.00873162225,   0.00963429082,    0.0106298644,    0.0117278378, 
	    0.012938669,    0.0142738689,    0.0157461073,    0.0173693337, 
	   0.0191588905,    0.0211316627,    0.0233062208,    0.0257029906, 
	   0.0283444282,    0.0312552229,    0.0344625004,    0.0379960723, 
	   0.0418886915,    0.0461763106,    0.0508984104,    0.0560983121, 
	   0.0618235618,    0.0681262836,    0.0750636458,    0.0826982856, 
	   0.0910988376,     0.100340441,     0.110505342,     0.121683538, 
	    0.133973435,     0.147482544,     0.162328333,     0.178639024, 
	    0.196554467,     0.216227174,     0.237823248,     0.261523515, 
	      0.2875247,     0.316040516,     0.347303122,       0.3815642, 
	    0.419096619,      0.46019569,     0.505180657,     0.554396272, 
	    0.608214319,     0.667035162,     0.731289208,     0.801438451, 
	    0.877978027,     0.961437345,      1.05238152,      1.15141225, 
	     1.25916874,      1.37632847,       1.5036068,      1.64175737, 
	     1.79157126,       1.9538753,      2.12953067,       2.3194294, 
	     2.52449131,       2.7456584,       2.9838903,       3.2401545, 
	     3.51541924,      3.81064296,      4.12676048,       4.4646697, 
	     4.82521486,        5.209167,      5.61720324,      6.04988337, 
	      6.5076232,      6.99066544,      7.49904966,      8.03257942, 
	     8.59078693,      9.17289543,      9.77778435,      10.4039488, 
	      11.049469,      11.7119675,      12.3885832,      13.0759459, 
	      13.770153,      14.4667568,      15.1607666,      15.8466597, 
	     16.5184078,      17.1695175,      17.7931042,      18.3819656, 
	     18.9287033,      19.4258461,      19.8660202,      20.2421207, 
	     20.5475311,      20.7763424,      20.9236012,      20.9855614, 
	      20.959938,      20.8461552,      20.6455784,      20.3617077, 
	     20.0003319,      19.5696163,      19.0801315,      18.5447674, 
	     17.9785786,      17.3984985,      16.8229542,      16.2713585, 
	     15.7635088,      15.3188791,      14.9558506,      14.6908951, 
	     14.5377407,      14.5065765,       14.603323,      14.8290148, 
	     15.1793613,      15.6444893,       16.208931,      16.8518715, 
	     17.5476608,      18.2665939,      18.9759541,      19.6412621, 
	     20.2276955,      20.7016373,      21.0322399,      21.1929798, 
	      21.163065,      20.9286613,      20.4838295,      19.8311119, 
	     18.9817276,      17.9553223,       16.779274,      15.4875593, 
	     14.1192303,      12.7165642,      11.3229771,      9.98083973, 
	     8.72929573,      7.60225439,      6.62667513,      5.82126474, 
	     5.19569063,      4.75036097,      4.47678614,       4.3584919, 
	     4.37240887,       4.4906168,      4.68230772,       4.9158144, 
	     5.16054201,      5.38867569,      5.57655668,      5.70564699, 
	     5.76306581,      5.74169493,      5.63989687,      5.46092463, 
	     5.21208477,      4.90376902,      4.54842377,      4.15954733, 
	       3.750772,      3.33508945,      2.92423534,       2.5282557, 
	     2.15524697,      1.81125629,      1.50031841,      1.22460055, 
	    0.984626889,     0.779551446,     0.607453465,     0.465633929, 
	    0.350893915,     0.259782642,     0.188806787,     0.134596705, 
	   0.0940291137,    0.0643085167,    0.0430115573,     0.028100213, 
	   0.0179103538,    0.0111223357,   0.00672012288,    0.0039446149, 
	  0.00224594655,   0.00123835925,  0.000660077203,  0.000339511345, 
	 0.000168188824,  8.00856578e-05,  3.65777341e-05,  1.59892934e-05, 
	 6.67409404e-06,  2.65374138e-06,  1.00260695e-06,  3.58971192e-07, 
	 1.21462506e-07,  3.87275527e-08,  1.16005303e-08,   3.2541867e-09, 
	 8.52084625e-10,   2.0754376e-10,  4.68574658e-11, 
/* Z=19 -> K */
	              0,   0.00160090986,   0.00176702545,   0.00195034686, 
	  0.00215265132,   0.00237589888,    0.0026222507,   0.00289409002, 
	  0.00319404528,    0.0035250138,   0.00389019004,   0.00429309579, 
	  0.00473761233,   0.00522801885,   0.00576903112,   0.00636584545, 
	  0.00702418759,   0.00775036775,   0.00855133589,   0.00943474751, 
	   0.0104090367,    0.0114834914,    0.0126683349,    0.0139748268, 
	   0.0154153621,    0.0170035791,    0.0187544916,    0.0206846129, 
	   0.0228121132,    0.0251569748,    0.0277411733,    0.0305888653, 
	   0.0337266028,    0.0371835679,    0.0409918129,    0.0451865494, 
	   0.0498064309,    0.0548938997,    0.0604955181,    0.0666623712, 
	   0.0734504759,    0.0809212402,    0.0891419575,     0.098186329, 
	    0.108135067,     0.119076468,     0.131107166,     0.144332737, 
	    0.158868626,     0.174840868,        0.192387,     0.211657077, 
	     0.23281464,     0.256037831,     0.281520486,     0.309473515, 
	    0.340125918,     0.373726457,     0.410544842,     0.450873256, 
	      0.4950279,     0.543350637,      0.59621042,     0.654005051, 
	    0.717162848,     0.786144137,     0.861442983,     0.943588436, 
	     1.03314626,      1.13072026,      1.23695266,      1.35252607, 
	       1.478163,      1.61462665,      1.76272035,       1.9232868, 
	     2.09720683,      2.28539681,      2.48880529,      2.70841026, 
	     2.94521189,      3.20022798,      3.47448397,      3.76900506, 
	     4.08480358,      4.42286634,       4.7841382,      5.16950464, 
	       5.579772,      6.01564407,      6.47769642,      6.96635008, 
	      7.4818387,       8.0241766,      8.59312439,      9.18815041, 
	     9.80839157,      10.4526167,      11.1191874,      11.8060131, 
	     12.5105228,      13.2296324,      13.9597073,      14.6965561, 
	     15.4354095,       16.170927,      16.8972073,      17.6078186, 
	     18.2958546,      18.9540005,      19.5746288,      20.1499157, 
	     20.6719933,      21.1331196,       21.525877,      21.8433952, 
	     22.0795956,      22.2294559,      22.2892685,      22.2569275, 
	     22.1321774,      21.9168549,      21.6150951,      21.2334843, 
	     20.7811451,      20.2697372,      19.7133751,      19.1284237, 
	     18.5331802,      17.9474468,      17.3919563,      16.8877125, 
	     16.4552155,      16.1136208,      15.8798552,      15.7677097, 
	     15.7869978,      15.9427748,      16.2347012,      16.6565952, 
	     17.1961937,      17.8351974,      18.5495872,      19.3102531, 
	     20.0839081,      20.8342819,      21.5235577,       22.113987, 
	     22.5696163,      22.8580494,      22.9521561,      22.8316135, 
	     22.4842205,      21.9068565,      21.1060371,      20.0979996, 
	     18.9082661,      17.5706997,      16.1260395,      14.6200294, 
	     13.1011467,      11.6181393,      10.2174416,      8.94069099, 
	     7.82247162,      6.88847637,      6.15420628,      5.62432194, 
	     5.29269505,       5.1431675,      5.15094566,      5.28452969, 
	     5.50801325,      5.78356409,      6.07389736,      6.34454298, 
	     6.56575584,      6.71394205,        6.772542,      6.73234653, 
	     6.59128046,      6.35373449,      6.02953434,      5.63268042, 
	      5.1799674,      4.68960667,      4.17994165,       3.6683383, 
	     3.17030001,      2.69883275,       2.2640636,      1.87310314, 
	      1.5301168,      1.23656869,     0.991590142,     0.792429626, 
	    0.634935737,     0.514036953,     0.424183816,     0.359732181, 
	    0.315251946,     0.285756409,     0.266854376,     0.254832804, 
	    0.246681541,     0.240074098,     0.233317003,     0.225279927, 
	    0.215315476,     0.203175858,     0.188930139,     0.172885045, 
	    0.155510679,     0.137371987,     0.119067319,     0.101175249, 
	   0.0842110366,    0.0685943887,    0.0546291657,    0.0424957424, 
	   0.0322547443,    0.0238607489,    0.0171831399,      0.01203107, 
	  0.00817923527,   0.00539160939,   0.00344093936,   0.00212279591, 
	  0.00126384653,  0.000724899641,  0.000399815792, 
/* Z=20 -> Ca */
	              0,   0.00189792144,   0.00209460221,   0.00231162622, 
	   0.0025510916,   0.00281531201,   0.00310683809,   0.00342848175, 
	  0.00378334359,    0.0041748411,   0.00460674195,   0.00508319773, 
	   0.0056087845,   0.00618854584,   0.00682803756,   0.00753338262, 
	  0.00831132568,   0.00916929543,    0.0101154773,    0.0111588798, 
	   0.0123094311,    0.0135780536,    0.0149767781,    0.0165188406, 
	   0.0182188135,    0.0200927258,    0.0221582167,    0.0244346838, 
	   0.0269434601,    0.0297080074,    0.0327541046,    0.0361100957, 
	   0.0398071148,    0.0438793637,    0.0483644083,    0.0533034876, 
	   0.0587418601,    0.0647291914,    0.0713199526,    0.0785738751, 
	   0.0865564197,    0.0953393206,     0.105001129,     0.115627863, 
	    0.127313614,     0.140161306,     0.154283434,     0.169802964, 
	    0.186854094,     0.205583319,     0.226150379,     0.248729363, 
	     0.27350992,      0.30069831,     0.330518931,     0.363215506, 
	     0.39905262,      0.43831709,     0.481319636,     0.528396547, 
	    0.579911172,     0.636255622,      0.69785279,       0.7651577, 
	    0.838659286,      0.91888231,      1.00638855,      1.10177886, 
	     1.20569384,      1.31881511,      1.44186676,      1.57561529, 
	     1.72086978,      1.87848198,      2.04934478,      2.23439121, 
	     2.43459129,      2.65095019,       2.8845017,       3.1363039, 
	     3.40743232,      3.69897056,      4.01200056,      4.34758949, 
	     4.70677757,      5.09055805,      5.49986219,      5.93553448, 
	     6.39831018,      6.88878679,      7.40739727,      7.95437479, 
	     8.52971745,      9.13315487,      9.76410484,      10.4216337, 
	     11.1044159,      11.8106918,       12.538229,      13.2842798, 
	     14.0455542,      14.8181849,      15.5977116,      16.3790703, 
	     17.1565914,      17.9240208,      18.6745453,      19.4008656, 
	     20.0952492,      20.7496548,      21.3558426,      21.9055481, 
	     22.3906517,      22.8034039,       23.136652,      23.3841152, 
	     23.5406475,      23.6025314,      23.5677586,      23.4363174, 
	     23.2104263,      22.8947582,      22.4965935,      22.0259018, 
	     21.4953403,      20.9201469,      20.3179073,      19.7082195, 
	     19.1121979,      18.5518684,      18.0494423,      17.6264725, 
	      17.302948,      17.0963326,      17.0206032,      17.0853424, 
	     17.2949009,      17.6477413,      18.1359768,      18.7451344, 
	     19.4542561,      20.2362728,      21.0587349,      21.8848686, 
	     22.6749153,      23.3877544,      23.9827061,       24.421463, 
	     24.6700363,      24.7006454,      24.4933968,      24.0376854, 
	     23.3331985,      22.3904362,      21.2306728,      19.8853416, 
	     18.3948097,      16.8065948,      15.1730738,      13.5488186, 
	     11.9876881,      10.5398579,      9.24900341,      8.14980793, 
	     7.26601219,      6.60914755,      6.17807627,      5.95937729, 
	     5.92856407,      6.05204439,      6.28965616,      6.59757996, 
	     6.93138027,      7.24894047,      7.51305008,      7.69347715, 
	     7.76837873,      7.72500181,      7.55967093,       7.2771349, 
	      6.8893795,      6.41405725,      5.87269688,      5.28884602, 
	     4.68629599,      4.08751631,      3.51237226,      2.97719741, 
	     2.49422908,      2.07140684,       1.7124933,      1.41746485, 
	     1.18310189,      1.00370657,     0.871871829,     0.779237151, 
	    0.717175186,     0.677370071,      0.65226388,     0.635364115, 
	    0.621417046,     0.606463492,     0.587798774,     0.563860238, 
	    0.534066916,     0.498632729,     0.458369374,     0.414493918, 
	    0.368450671,     0.321754307,     0.275858998,     0.232056737, 
	    0.191405147,     0.154685274,     0.122387007,    0.0947189406, 
	   0.0716383681,    0.0528957993,    0.0380881317,    0.0267143175, 
	   0.0182281937,     0.012084079,   0.00777220773,   0.00484274188, 
	  0.00291859009,    0.0016985205,  0.000952854462,  0.000514330168, 
	 0.000266609801,   0.00013244756,  6.29237111e-05, 
/* Z=21 -> Sc */
	              0,   0.00223217905,   0.00246318663,   0.00271805376, 
	  0.00299923681,   0.00330944429,   0.00365166203,   0.00402918132, 
	  0.00444562966,   0.00490500545,   0.00541171478,   0.00597061403, 
	  0.00658705365,   0.00726692798,   0.00801673159,   0.00884361845, 
	    0.009755468,    0.0107609574,     0.011869642,    0.0130920429, 
	    0.014439743,    0.0159254912,     0.017563317,    0.0193686616, 
	   0.0213585123,    0.0235515516,    0.0259683337,    0.0286314543, 
	   0.0315657593,    0.0347985551,    0.0383598469,    0.0422826037, 
	   0.0466030352,    0.0513609014,    0.0565998517,    0.0623677894, 
	   0.0687172636,    0.0757059082,    0.0833969191,    0.0918595418, 
	    0.101169661,     0.111410342,     0.122672543,     0.135055751, 
	    0.148668811,     0.163630635,     0.180071175,     0.198132306, 
	    0.217968792,     0.239749461,     0.263658255,     0.289895445, 
	    0.318679065,     0.350245982,     0.384853631,     0.422781378, 
	    0.464332104,     0.509833872,      0.55964154,     0.614138722, 
	    0.673739254,     0.738889277,     0.810068905,     0.887794137, 
	    0.972618282,      1.06513417,      1.16597509,       1.2758168, 
	     1.39537799,      1.52542186,      1.66675651,      1.82023454, 
	     1.98675346,      2.16725445,      2.36272001,       2.5741725, 
	     2.80266976,      3.04930091,      3.31518006,        3.601439, 
	     3.90921807,      4.23965597,      4.59387541,      4.97296906, 
	     5.37798309,      5.80989456,       6.2695899,      6.75784159, 
	     7.27527571,       7.8223443,      8.39928913,      9.00610447, 
	     9.64250088,      10.3078613,      11.0011988,      11.7211132, 
	     12.4657469,      13.2327461,      14.0192146,      14.8216848, 
	     15.6360865,      16.4577198,      17.2812576,      18.1007328, 
	     18.9095669,      19.7006016,      20.4661598,      21.1981258, 
	     21.8880501,      22.5272999,      23.1072006,      23.6192608, 
	     24.0553703,      24.4080734,      24.6708298,       24.838316, 
	     24.9067192,      24.8740482,      24.7404118,      24.5082932, 
	     24.1827679,      23.7716579,      23.2856312,      22.7381802, 
	     22.1455097,      21.5262928,      20.9012814,      20.2928085, 
	     19.7241135,      19.2185688,      18.7987881,      18.4856453, 
	     18.2972507,      18.2479305,      18.3472519,      18.5991459, 
	     19.0012112,      19.5442238,      20.2119255,      20.9811325, 
	     21.8221855,      22.6997585,      23.5740185,      24.4021358, 
	     25.1400337,      25.7443962,      26.1747646,      26.3956718, 
	     26.3786793,      26.1041908,      25.5629387,      24.7570171, 
	     23.7003574,      22.4185886,      20.9482365,      19.3352547, 
	     17.6329155,      15.8991737,      14.1936321,      12.5742521, 
	     11.0940733,      9.79811859,      8.72075462,      7.88368988, 
	     7.29480267,        6.947896,      6.82342768,      6.89016294, 
	     7.10762119,      7.42912102,        7.805161,      8.18685627, 
	     8.52915382,      8.79357338,      8.95027351,      8.97933769, 
	     8.87121391,      8.62636948,       8.2542429,      7.77165556, 
	     7.20086622,      6.56746149,      5.89825583,      5.21937513, 
	     4.55463505,      3.92430329,      3.34428549,       2.8257339, 
	     2.37505245,      1.99423218,      1.68144488,      1.43180525, 
	     1.23821485,      1.09221101,     0.984749019,      0.90687418, 
	     0.85025084,     0.807538986,     0.772621512,     0.740700424, 
	    0.708284676,     0.673096716,     0.633924007,       0.5904392, 
	    0.543006599,      0.49249053,     0.440075576,     0.387104928, 
	    0.334942192,     0.284859151,     0.237950429,     0.195075914, 
	    0.156829625,     0.123533197,     0.095250845,    0.0718215555, 
	   0.0529035032,     0.038025178,    0.0266378131,    0.0181643143, 
	   0.0120407473,   0.00774801103,   0.00483263843,   0.00291712931, 
	   0.0017013289,  0.000957036391,  0.000518302142,  0.000269722776, 
	 0.000134603106,  6.42795203e-05,   2.9308947e-05, 
/* Z=22 -> Ti */
	              0,   0.00260778167,   0.00287727988,   0.00317457109, 
	  0.00350251189,   0.00386425178,   0.00426326133,   0.00470336573, 
	   0.0051887813,   0.00572415255,   0.00631459849,   0.00696575874, 
	   0.0076838457,   0.00847570319,   0.00934886746,    0.0103116408, 
	   0.0113731651,    0.0125435032,     0.013833737,    0.0152560631, 
	   0.0168239046,    0.0185520351,    0.0204567071,    0.0225558039, 
	   0.0248689894,    0.0274178907,    0.0302262884,    0.0333203189, 
	   0.0367287174,    0.0404830426,    0.0446179733,    0.0491715893, 
	   0.0541857108,     0.059706226,    0.0657835007,    0.0724727809, 
	   0.0798346549,    0.0879355371,    0.0968482122,     0.106652409, 
	    0.117435433,     0.129292816,     0.142329127,     0.156658635, 
	    0.172406316,     0.189708591,      0.20871447,     0.229586497, 
	    0.252501905,     0.277653813,     0.305252522,     0.335526824, 
	    0.368725389,     0.405118406,     0.444999069,     0.488685191, 
	    0.536521018,      0.58887893,     0.646161318,     0.708802462, 
	    0.777270317,     0.852068543,     0.933738291,      1.02286005, 
	     1.12005556,       1.2259891,      1.34136915,       1.4669497, 
	     1.60353112,      1.75196052,      1.91313255,       2.0879879, 
	     2.27751398,      2.48274136,      2.70474219,      2.94462562, 
	      3.2035346,      3.48263741,      3.78312206,      4.10618544, 
	       4.453022,      4.82481098,       5.2227006,      5.64778852, 
	     6.10110378,      6.58357906,      7.09602833,      7.63911295, 
	      8.2133131,      8.81888771,      9.45584202,      10.1238794, 
	     10.8223629,      11.5502701,      12.3061457,      13.0880594, 
	     13.8935547,      14.7196159,      15.5626249,      16.4183331, 
	     17.2818394,      18.1475735,      19.0093079,      19.8601665, 
	     20.6926689,      21.4987907,      22.2700424,      22.9975948, 
	     23.6724091,      24.2854195,      24.8277359,      25.2908764, 
	     25.6670361,      25.9493752,      26.1323242,      26.2119026, 
	     26.1860352,      26.0548553,      25.8209877,      25.4897766, 
	     25.0694618,      24.5712605,      24.0093689,      23.4008312, 
	     22.7652817,      22.1245441,      21.5020905,      20.9223404, 
	     20.4098377,      19.9882946,      19.6795712,       19.502573, 
	     19.4721699,      19.5981579,      19.8843346,      20.3277493, 
	       20.91819,      21.6379642,      22.4620209,      23.3584328, 
	     24.2893009,      25.2119904,      26.0807781,      26.8487606, 
	     27.4700108,      27.9018669,      28.1072369,      28.0567932, 
	     27.7309322,      27.1213398,      26.2320709,      25.0800095, 
	     23.6946297,      22.1170425,      20.3982754,      18.5968933, 
	     16.7760162,      14.9999237,        13.33043,      11.8232737, 
	     10.5247898,      9.46910095,      8.67608833,      8.15029526, 
	     7.88089561,      7.84274626,       7.9984436,      8.30122948, 
	     8.69849586,      9.13559628,      9.55963326,      9.92291451, 
	     10.1858158,      10.3188295,      10.3037157,      10.1337042, 
	     9.81283665,      9.35457993,      8.77988625,      8.11493683, 
	     7.38877487,      6.63102198,      5.86985826,      5.13038445, 
	     4.43342829,      3.79483891,      3.22523808,      2.73018169, 
	     2.31065226,      1.96379077,      1.68376637,      1.46269345, 
	     1.29151475,      1.16078866,      1.06133795,     0.984740257, 
	    0.923657775,     0.872019708,     0.825078726,      0.77937007, 
	    0.732600987,      0.68349576,     0.631618083,     0.577186644, 
	    0.520896018,     0.463749707,     0.406911761,     0.351578653, 
	    0.298874706,     0.249771237,      0.20503059,     0.165173918, 
	     0.13047199,     0.100955911,    0.0764447376,    0.0565851778, 
	   0.0408986397,    0.0288303923,    0.0197962895,     0.013223161, 
	  0.00858027115,   0.00540063996,   0.00329229725,   0.00194070255, 
	  0.00110429549,  0.000605480571,  0.000319290615,  0.000161615724, 
	 7.83590658e-05,  3.63122963e-05,  1.60464515e-05, 
/* Z=23 -> V */
	              0,   0.00302829151,   0.00334078679,    0.0036854588, 
	  0.00406560861,   0.00448487559,   0.00494726934,   0.00545720849, 
	  0.00601956155,   0.00663969247,   0.00732350908,   0.00807752088, 
	  0.00890889671,   0.00982553046,     0.010836117,    0.0119502284, 
	   0.0131784016,     0.014532241,    0.0160245113,    0.0176692642, 
	   0.0194819625,    0.0214796141,    0.0236809272,    0.0261064824, 
	   0.0287789013,    0.0317230597,    0.0349663012,    0.0385386646, 
	   0.0424731597,    0.0468060486,    0.0515771471,     0.056830164, 
	     0.06261307,    0.0689785182,    0.0759842396,     0.083693549, 
	   0.0921758413,     0.101507172,     0.111770831,     0.123058014, 
	    0.135468528,     0.149111569,     0.164106503,      0.18058382, 
	    0.198685989,     0.218568593,     0.240401343,     0.264369279, 
	    0.290673971,     0.319534987,     0.351191103,     0.385901928, 
	     0.42394945,     0.465639591,      0.51130414,     0.561302245, 
	    0.616022587,     0.675885022,     0.741342723,     0.812884033, 
	    0.891034544,      0.97635895,      1.06946301,      1.17099547, 
	     1.28164959,      1.40216494,      1.53332818,      1.67597485, 
	     1.83098948,      1.99930584,      2.18190622,      2.37982178, 
	     2.59412885,      2.82594728,      3.07643723,      3.34679222, 
	     3.63823485,      3.95200729,      4.28936243,      4.65155172, 
	     5.03981066,       5.4553442,      5.89930582,      6.37277699, 
	     6.87674284,      7.41206312,      7.97944307,      8.57939816, 
	     9.21221733,      9.87792492,      10.5762329,      11.3065023, 
	     12.0676918,      12.8583136,      13.6763811,      14.5193663, 
	     15.3841572,      16.2670135,      17.1635418,      18.0686626, 
	     18.9766026,      19.8808956,      20.7744045,      21.6493511, 
	     22.4973888,      23.3096771,      24.0770149,      24.7899723, 
	     25.4390869,      26.0150566,      26.5090065,      26.9127483, 
	     27.2190914,      27.4221535,       27.517704,      27.5034904, 
	     27.3795605,      27.1485558,      26.8159618,      26.3902836, 
	     25.8831558,      25.3093338,      24.6865635,      24.0353241, 
	     23.3784142,      22.7403774,        22.14678,      21.6233425, 
	     21.1949444,      20.8845329,      20.7119732,      20.6928978, 
	     20.8375988,      21.1500511,      21.6271038,      22.2579231, 
	     23.0237637,      23.8980713,      24.8470116,      25.8304043, 
	     26.8030624,      27.7165375,      28.5211658,      29.1683979, 
	     29.6132412,      29.8167515,      29.7483845,      29.3881073, 
	     28.7280693,      27.7737484,      26.5444031,      25.0727825, 
	     23.4040222,      21.5937309,      19.7053356,      17.8067989, 
	     15.9668694,      14.2511311,      12.7180834,      11.4155521, 
	     10.3777332,      9.62308979,      9.15332603,      8.95353603, 
	       8.993536,      9.23029613,      9.61126423,      10.0783033, 
	     10.5718966,      11.0352736,      11.4180984,      11.6794453, 
	      11.789856,      11.7323732,      11.5025368,      11.1074638, 
	      10.564147,      9.89722347,      9.13643742,      8.31404495, 
	     7.46236992,      6.61169434,      5.78858757,      5.01475143, 
	     4.30638981,      3.67406535,      3.12298203,      2.65359211, 
	     2.26243043,      1.94306135,      1.68705046,      1.48487425, 
	     1.32671237,      1.20308673,      1.10533011,      1.02589345, 
	    0.958505213,     0.898211241,     0.841324329,     0.785310566, 
	    0.728638411,     0.670608819,     0.611181259,     0.550805449, 
	    0.490264446,     0.430533499,     0.372656882,     0.317644536, 
	    0.266390145,     0.219611004,     0.177810326,     0.141260996, 
	     0.11000926,    0.0838954151,    0.0625875592,    0.0456240922, 
	   0.0324599445,    0.0225119591,    0.0151994573,   0.00997701101, 
	  0.00635782769,   0.00392734911,   0.00234793499,   0.00135627261, 
	 0.000755655812,  0.000405340776,  0.000208928934,  0.000103271093, 
	 4.88467667e-05,   2.2059623e-05,  9.48951947e-06, 
/* Z=24 -> Cr */
	              0,   0.00349577097,   0.00385595392,     0.004253163, 
	  0.00469119055,   0.00517421588,   0.00570684299,   0.00629414432, 
	  0.00694170734,   0.00765568763,   0.00844286382,   0.00931070186, 
	   0.0102674216,    0.0113220774,    0.0124846334,    0.0137660578, 
	    0.015178428,    0.0167350303,    0.0184504874,    0.0203408841, 
	   0.0224239212,    0.0247190576,    0.0272476971,    0.0300333723, 
	   0.0331019498,    0.0364818648,    0.0402043462,     0.044303719, 
	   0.0488176718,    0.0537875853,    0.0592589006,    0.0652814731, 
	   0.0719099939,    0.0792044699,    0.0872306749,    0.0960607007, 
	    0.105773553,     0.116455749,     0.128202006,     0.141115978, 
	    0.155311063,     0.170911208,     0.188051894,     0.206881076, 
	    0.227560237,     0.250265598,     0.275189251,     0.302540421, 
	    0.332547039,     0.365456909,     0.401539564,     0.441087633, 
	     0.48441869,     0.531876981,     0.583835483,      0.64069742, 
	    0.702898622,     0.770909488,     0.845236957,      0.92642653, 
	     1.01506436,      1.11177957,      1.21724558,      1.33218229, 
	     1.45735753,      1.59358883,      1.74174368,      1.90274143, 
	     2.07755232,      2.26719761,      2.47274923,      2.69532681, 
	     2.93609619,      3.19626379,      3.47707391,      3.77979994, 
	     4.10573769,       4.4561944,      4.83247709,      5.23587942, 
	     5.66766214,       6.1290369,      6.62114096,      7.14501476, 
	     7.70157051,      8.29156208,      8.91555119,      9.57386589, 
	     10.2665625,      10.9933758,      11.7536821,      12.5464401, 
	     13.3701496,      14.2227964,      15.1018047,      16.0039978, 
	     16.9255428,      17.8619328,      18.8079472,      19.7576466, 
	     20.7043629,      21.6407261,      22.5587025,       23.449646, 
	     24.3043976,      25.1134033,      25.8668575,      26.5548954, 
	     27.1678104,      27.6963062,      28.1317787,      28.4666309, 
	     28.6946087,      28.8111439,      28.8136997,      28.7021179, 
	     28.4789181,      28.1495667,      27.7226772,      27.2101212, 
	     26.6270351,      25.9916973,      25.3252754,      24.6513901, 
	     23.9955502,      23.3843918,      22.8447876,      22.4028034, 
	     22.0825672,      21.9050579,      21.8868961,      22.0391712, 
	     22.3663979,      22.8656502,      23.5259609,      24.3280334, 
	     25.2443562,      26.2397156,      27.2721577,      28.2943993, 
	     29.2556458,      30.1037617,      30.7877197,      31.2602119, 
	     31.4802856,      31.4158611,      31.0459824,      30.3626003, 
	     29.3717995,      28.0942726,      26.5650101,      24.8320904, 
	     22.9546223,      20.9998837,      19.0397873,       17.146862, 
	     15.3900204,      13.8303881,      12.5175056,       11.486227, 
	     10.7545624,      10.3226824,      10.1731901,      10.2726593, 
	     10.5743237,      11.0216951,      11.5527906,      12.1046104, 
	     12.6174631,      13.0387983,      13.3262367,      13.4496088, 
	     13.3918886,       13.149066,      12.7290545,      12.1498365, 
	     11.4370832,      10.6215115,      9.73622513,       8.8142643, 
	     7.88650894,      6.98007774,      6.11724663,      5.31489372, 
	     4.58442593,       3.9321022,      3.35965347,      2.86509776, 
	     2.44363904,      2.08856583,      1.79207587,      1.54597723, 
	     1.34223986,      1.17339075,      1.03276086,     0.914607465, 
	    0.814136565,     0.727456212,     0.651486516,     0.583849609, 
	    0.522755265,     0.466893524,     0.415339351,     0.367471337, 
	    0.322902918,     0.281424493,     0.242953405,     0.207490489, 
	    0.175081566,      0.14578405,     0.119638793,    0.0966482013, 
	   0.0767613202,    0.0598662794,    0.0457901582,    0.0343054123, 
	   0.0251416303,    0.0180007499,    0.0125737227,   0.00855665654, 
	  0.00566475047,   0.00364283333,   0.00227194838,   0.00137198356, 
	 0.000800846552,  0.000451050902,  0.000244662486,  0.000127563457, 
	  6.3798856e-05,  3.05418507e-05,  1.39635513e-05, 
/* Z=25 -> Mn */
	              0,   0.00401948253,   0.00443296554,   0.00488888146, 
	  0.00539156888,   0.00594580779,   0.00655686157,   0.00723053003, 
	  0.00797319971,   0.00879190397,   0.00969438814,    0.0106891794, 
	   0.0117856665,    0.0129941842,    0.0143261068,    0.0157939568, 
	   0.0174115058,    0.0191939119,    0.0211578496,    0.0233216584, 
	    0.025705507,    0.0283315722,    0.0312242359,    0.0344102904, 
	   0.0379191898,    0.0417832844,    0.0460381135,    0.0507227108, 
	   0.0558799282,    0.0615568012,    0.0678049475,     0.074680984, 
	   0.0822470114,    0.0905711129,    0.0997278988,     0.109799117, 
	    0.120874301,      0.13305144,     0.146437794,     0.161150664, 
	    0.177318275,     0.195080757,     0.214591175,     0.236016557, 
	    0.259539157,     0.285357624,     0.313688397,     0.344767153, 
	    0.378850222,     0.416216284,     0.457168013,      0.50203383, 
	    0.551169813,     0.604961574,     0.663826287,     0.728214741, 
	    0.798613489,     0.875546813,     0.959579229,      1.05131733, 
	     1.15141165,      1.26055956,      1.37950599,      1.50904608, 
	      1.6500262,      1.80334508,      1.96995485,      2.15086126, 
	     2.34712291,      2.55985165,      2.79020929,      3.03940535, 
	     3.30869365,      3.59936714,      3.91275048,      4.25019312, 
	     4.61305666,      5.00270605,      5.42049217,      5.86773682, 
	     6.34571028,      6.85561085,      7.39853811,      7.97546387, 
	     8.58719826,      9.23435593,      9.91731453,      10.6361732, 
	      11.390707,       12.180316,      13.0039797,      13.8602028, 
	     14.7469616,      15.6616592,      16.6010666,      17.5612946, 
	     18.5377407,      19.5250721,      20.5172024,      21.5072918, 
	     22.4877644,      23.4503422,      24.3861065,      25.2855873, 
	     26.1388855,      26.9358253,      27.6661358,      28.3196831, 
	     28.8867264,      29.3582115,      29.7260895,      29.9836712, 
	     30.1259766,      30.1501122,      30.0556126,      29.8447781, 
	     29.5229588,      29.0987625,      28.5841999,      27.9946938, 
	     27.3489799,      26.6688442,      25.9787102,      25.3050404, 
	     24.6755791,      24.1184254,      23.6609535,      23.3286266, 
	     23.1437397,      23.1241322,      23.2819462,      23.6225014, 
	     24.1433563,       24.833643,      25.6737213,      26.6352444, 
	     27.6816673,      28.7692242,       29.848381,      30.8657284, 
	     31.7662811,      32.4960632,      33.0048943,      33.2492256, 
	     33.1948433,      32.8193054,      32.1139145,      31.0850353, 
	     29.7546864,      28.1601906,      26.3529339,      24.3961353, 
	     22.3617649,      20.3267078,      18.3684025,      16.5602283, 
	      14.966958,      13.6406317,      12.6171789,      11.9141016, 
	     11.5294342,      11.4421043,      11.6136951,      11.9914894, 
	     12.5125141,      13.1082678,      13.7096863,      14.2519331, 
	     14.6786098,      14.9450502,       15.020484,      14.8889732, 
	     14.5491371,      14.0128069,      13.3028402,      12.4503632, 
	     11.4917297,      10.4655018,      9.40966225,      8.35927773, 
	     7.34471416,      6.39046288,       5.5145731,      4.72863245, 
	     4.03820229,      3.44358754,      2.94082403,      2.52275705, 
	     2.18011284,      1.90247881,      1.67913735,      1.49972415, 
	     1.35470271,      1.23566759,      1.13550031,      1.04840708, 
	      0.9698717,      0.89654994,     0.826131821,     0.757187486, 
	    0.689009488,     0.621459424,     0.554822564,     0.489673823, 
	    0.426756322,     0.366874665,     0.310804486,     0.259219676, 
	    0.212639153,     0.171393409,     0.135610744,     0.105221599, 
	   0.0799785256,    0.0594878756,    0.0432493016,    0.0306981467, 
	   0.0212466568,    0.0143202916,    0.0093865525,   0.00597495493, 
	  0.00368797616,   0.00220386451,   0.00127294369,  0.000709429325, 
	 0.000380798127,  0.000196488094,  9.72659982e-05,  4.60948941e-05, 
	 2.08664333e-05,  9.00192572e-06,  3.69191321e-06, 
/* Z=26 -> Fe */
	              0,   0.00459854677,   0.00507081486,    0.0055914633, 
	  0.00616542948,   0.00679815048,   0.00749561517,   0.00826441869, 
	  0.00911182258,    0.0100458227,    0.0110752219,    0.0122097088, 
	   0.0134599488,    0.0148376813,    0.0163558219,    0.0180285797, 
	   0.0198715925,    0.0219020564,    0.0241388846,    0.0266028736, 
	    0.029316891,    0.0323060714,    0.0355980359,     0.039223142, 
	   0.0432147346,    0.0476094373,    0.0524474718,    0.0577729903, 
	   0.0636344478,    0.0700850263,    0.0771830454,    0.0849924684, 
	   0.0935834125,     0.103032723,     0.113424577,      0.12485116, 
	    0.137413383,     0.151221663,     0.166396767,     0.183070704, 
	    0.201387748,     0.221505418,      0.24359569,     0.267846137, 
	    0.294461221,     0.323663771,     0.355696321,     0.390822709, 
	    0.429329813,     0.471529156,     0.517758846,     0.568385303, 
	    0.623805642,     0.684449315,     0.750780404,     0.823299885, 
	    0.902547777,     0.989105403,      1.08359754,      1.18669474, 
	     1.29911506,      1.42162669,      1.55504918,      1.70025527, 
	     1.85817218,      2.02978206,       2.2161231,      2.41828895, 
	     2.63742781,      2.87474084,       3.1314795,      3.40894127, 
	     3.70846581,      4.03142691,       4.3792243,      4.75327492, 
	     5.15499926,      5.58580542,      6.04707527,      6.54014158, 
	     7.06626415,      7.62660599,       8.2222023,      8.85392666, 
	     9.52245426,      10.2282219,       10.971386,        11.75177, 
	      12.568821,      13.4215536,      14.3085012,      15.2276554, 
	     16.1764183,      17.1515522,       18.149128,      19.1644917, 
	     20.1922283,      21.2261467,      22.2592716,      23.2838554, 
	     24.2914143,      25.2727928,      26.2182388,      27.1175385, 
	     27.9601536,       28.735424,      29.4327888,      30.0420456, 
	     30.5536613,      30.9590969,      31.2511578,      31.4243774, 
	      31.475399,      31.4033337,      31.2101231,      30.9008369, 
	     30.4839134,      29.9713135,      29.3785648,      28.7246666, 
	     28.0318527,      27.3251762,      26.6319313,      25.9808693, 
	     25.4012604,      24.9217949,      24.5693512,      24.3676796, 
	     24.3360615,      24.4879894,      24.8299713,      25.3605042, 
	     26.0693245,      26.9369984,      27.9349327,       29.025835, 
	     30.1646805,      31.3001823,      32.3767395,      33.3368149, 
	     34.1236382,      34.6841583,      34.9720192,       34.950489, 
	     34.5950623,      33.8956108,      32.8578568,      31.5040646, 
	     29.8727589,      28.0174465,      26.0043221,      23.9090176, 
	     21.8125362,      19.7966232,      17.9388237,       16.307621, 
	     14.9579964,      13.9277925,      13.2352057,      12.8776522, 
	     12.8321428,      13.0571766,      13.4960003,       14.080965, 
	     14.7386084,      15.3949957,      15.9808712,      16.4361629, 
	     16.7135105,      16.7805672,      16.6209583,      16.2339516, 
	     15.6329765,       14.843236,      13.8987226,      12.8389406, 
	     11.7056589,      10.5399389,       9.3796339,        8.257514, 
	     7.20002842,      6.22672415,      5.35024834,       4.5768218, 
	      3.9070704,      3.33707309,      2.85950398,      2.46476007, 
	     2.14198899,      1.87996125,      1.66775692,      1.49526179, 
	     1.35348618,      1.23473144,       1.1326369,      1.04213989, 
	    0.959377646,     0.881555974,     0.806801796,     0.734011948, 
	    0.662704527,     0.592877388,     0.524875998,     0.459271461, 
	    0.396751195,     0.338023424,     0.283737689,     0.234423116, 
	    0.190445587,     0.151984483,     0.119027764,    0.0913837254, 
	   0.0687063783,    0.0505302213,    0.0363102779,    0.0254629273, 
	   0.0174036287,    0.0115784509,   0.00748754106,   0.00469976198, 
	  0.00285890047,   0.00168273016,  0.000956738601,  0.000524529256, 
	 0.000276782812,  0.000140299555,  6.81761885e-05,  3.16909282e-05, 
	 1.40598959e-05,  5.93941468e-06,  2.38307257e-06, 
/* Z=27 -> Co */
	              0,   0.00523932744,   0.00577647844,   0.00636855653, 
	  0.00702115521,   0.00774043472,    0.0085331779,    0.0094068516, 
	   0.0103696771,     0.011430705,    0.0125998948,    0.0138882091, 
	   0.0153077114,     0.016871674,     0.018594699,    0.0204928499, 
	   0.0225837883,    0.0248869453,    0.0274236798,    0.0302174706, 
	   0.0332941264,    0.0366820134,    0.0404122844,    0.0445191748, 
	   0.0490402728,    0.0540168583,    0.0594942421,    0.0655221418, 
	   0.0721551254,    0.0794530213,     0.087481454,    0.0963123515, 
	    0.106024519,       0.1167043,     0.128446251,     0.141353846, 
	    0.155540317,      0.17112948,     0.188256741,     0.207069978, 
	    0.227730751,     0.250415325,      0.27531603,     0.302642554, 
	    0.332623273,     0.365506858,     0.401563883,     0.441088408, 
	    0.484399885,     0.531844974,     0.583799541,     0.640670657, 
	     0.70289892,     0.770960271,     0.845368683,     0.926678181, 
	     1.01548529,      1.11243129,      1.21820426,      1.33354199, 
	     1.45923305,      1.59611952,      1.74509811,      1.90712214, 
	     2.08320165,        2.274405,      2.48185778,      2.70674253, 
	     2.95029759,       3.2138133,      3.49862957,      3.80612993, 
	     4.13773537,      4.49489641,       4.8790822,      5.29176712, 
	     5.73441982,      6.20848083,      6.71534538,      7.25634003, 
	      7.8326931,        8.445508,      9.09572697,      9.78409576, 
	      10.511117,      11.2770147,      12.0816774,      12.9246101, 
	     13.8048849,      14.7210789,      15.6712227,      16.6527462, 
	     17.6624203,      18.6963158,      19.7497501,      20.8172684, 
	     21.8925991,      22.9686661,      24.0375919,      25.0907173, 
	     26.1186771,      27.1114731,      28.0585899,      28.9491558, 
	     29.7721214,      30.5164986,      31.1716194,      31.7274437, 
	     32.1748962,      32.5062332,      32.7154274,      32.7985611, 
	     32.7542229,      32.5838661,      32.2921486,      31.8871803, 
	     31.3807201,       30.788229,      30.1288128,      29.4249973, 
	     28.7023411,      27.9888458,      27.3142014,      26.7088203, 
	     26.2027302,       25.824297,       25.598875,      25.5473976, 
	     25.6849842,      26.0196686,      26.5512714,      27.2705727, 
	     28.1588058,      29.1875896,      30.3193264,      31.5081367, 
	     32.7013092,      33.8412933,      34.8681068,      35.7221756, 
	     36.3473969,      36.6943092,      36.7231865,      36.4068489, 
	     35.7330284,      34.7060204,      33.3475151,      31.6964207, 
	     29.8076344,      27.7496986,      25.6014423,      23.4477348, 
	     21.3745918,      19.4639416,      17.7884159,      16.4065838, 
	     15.3590078,      14.6654844,      14.3237514,       14.309803, 
	     14.5798311,      15.0736513,      15.7193232,      16.4385548, 
	     17.1524277,       17.786911,      18.2777271,      18.5741501, 
	     18.6415024,      18.4622211,      18.0355301,      17.3758659, 
	     16.5103321,      15.4755087,      14.3139448,      13.0706749, 
	     11.7900248,      10.5129509,      9.27500343,      8.10501194, 
	     7.02444744,      6.04741526,      5.18114948,       4.4268856, 
	     3.78096581,      3.23604059,      2.78225827,      2.40834951, 
	     2.10254931,      1.85332477,      1.64990497,      1.48262453, 
	     1.34310997,      1.22434175,      1.12062478,      1.02750015, 
	    0.941620648,     0.860610127,      0.78291738,     0.707672119, 
	    0.634546578,     0.563624263,     0.495277762,     0.430056244, 
	     0.36858508,      0.31147939,     0.259273469,     0.212368339, 
	    0.170998216,     0.135215521,     0.104893707,    0.0797449127, 
	   0.0593491681,    0.0431910753,    0.0306994505,    0.0212858766, 
	   0.0143787032,   0.00945000257,   0.00603416376,   0.00373795768, 
	  0.00224290695,   0.00130149303,  0.000729096122,  0.000393605675, 
	 0.000204386917,  0.000101882673,  4.86518293e-05,  2.22075323e-05, 
	  9.6673175e-06,  4.00376075e-06,  1.57357044e-06, 
/* Z=28 -> Ni */
	              0,   0.00594680198,   0.00655539893,    0.0072261109, 
	  0.00796525087,   0.00877976976,   0.00967731792,    0.0106663182, 
	   0.0117560402,    0.0129566845,    0.0142794782,    0.0157367717, 
	   0.0173421577,    0.0191105809,    0.0210584812,    0.0232039429, 
	   0.0255668443,    0.0281690415,    0.0310345571,    0.0341897942, 
	   0.0376637653,    0.0414883383,    0.0456985272,     0.050332766, 
	   0.0554332621,     0.061046347,    0.0672228485,     0.074018538, 
	   0.0814945921,    0.0897180587,    0.0987624601,     0.108708329, 
	    0.119643912,     0.131665796,     0.144879743,     0.159401417, 
	    0.175357327,     0.192885801,     0.212137923,     0.233278707, 
	    0.256488264,      0.28196305,     0.309917241,     0.340584159, 
	    0.374217868,     0.411094695,     0.451515198,     0.495805651, 
	    0.544320285,     0.597443223,     0.655590415,     0.719212115, 
	    0.788794935,     0.864864171,     0.947986484,      1.03877175, 
	     1.13787603,      1.24600315,      1.36390829,      1.49239874, 
	      1.6323365,      1.78464055,      1.95028734,      2.13031316, 
	     2.32581306,      2.53794289,      2.76791692,      3.01700783, 
	     3.28654289,       3.5779016,      3.89250994,      4.23183393, 
	     4.59737253,      4.99064541,      5.41318274,      5.86650896, 
	     6.35212708,      6.87149572,      7.42600727,       8.0169611, 
	     8.64553165,      9.31273556,       10.019393,      10.7660866, 
	     11.5531158,      12.3804445,       13.247654,      14.1538858, 
	      15.097784,      16.0774384,      17.0903263,      18.1332626, 
	     19.2023373,      20.2928791,      21.3994122,      22.5156345, 
	     23.6343994,      24.7477283,      25.8468285,      26.9221592, 
	     27.9634972,      28.9600639,      29.9006634,      30.7738838, 
	      31.568306,      32.2727928,      32.8767776,      33.3706169, 
	     33.7459717,      33.9961777,      34.1166878,      34.1054573, 
	     33.9633408,      33.6944351,       33.306385,      32.8105774, 
	     32.2222404,      31.5604248,      30.8477783,      30.1102142, 
	     29.3763237,      28.6766357,      28.0426559,      27.5057392, 
	     27.0958042,       26.839922,      26.7608662,      26.8756466, 
	     27.1941376,      27.7178802,      28.4391289,      29.3402481, 
	     30.3935509,      31.5616112,      32.7981262,      34.0493393, 
	      35.256012,      36.3559036,      37.2866516,      37.9889832, 
	     38.4100304,      38.5066185,      38.2483139,       37.619976, 
	      36.623661,      35.2796249,      33.6263008,      31.7191753, 
	     29.6284466,       27.435627,      25.2291412,      23.0992146, 
	     21.1323471,      19.4057713,      17.9823151,      16.9061203, 
	     16.1995544,      15.8616972,      15.8684988,      16.1747093, 
	     16.7173824,      17.4207115,      18.2017422,      18.9764481, 
	      19.665659,      20.2002945,      20.5255299,      20.6035347, 
	     20.4147339,      19.9575291,      19.2467117,      18.3107891, 
	     17.1886139,      15.9256449,      14.5702276,      13.1701632, 
	     11.7698154,      10.4079189,      9.11612225,      7.91826677, 
	       6.830338,      5.86095715,      5.01227951,        4.281147, 
	     3.66035199,      3.13989115,      2.70811772,      2.35272336, 
	      2.0615201,      1.82301402,      1.62678516,      1.46370161, 
	     1.32600152,      1.20727861,      1.10240495,      1.00741386, 
	    0.919366837,     0.836211026,      0.75663954,     0.679955602, 
	    0.605942369,     0.534740031,     0.466730863,     0.402433902, 
	    0.342411727,     0.287191451,     0.237202317,     0.192731574, 
	     0.15389888,     0.120649122,    0.0927614346,    0.0698717237, 
	   0.0515046604,    0.0371112078,    0.0261071101,    0.0179088134, 
	    0.011963645,   0.00777243217,   0.00490372907,   0.00300000072, 
	  0.00177688233,   0.00101723499,   0.00056189613,  0.000298931554, 
	 0.000152875145,   7.5002863e-05,  3.52270581e-05,  1.58040202e-05, 
	 6.75674983e-06,  2.74613785e-06,  1.05828644e-06, 
/* Z=29 -> Cu */
	              0,   0.00672325259,   0.00741004106,   0.00816678908, 
	  0.00900059007,   0.00991925504,    0.0109313782,      0.01204642, 
	    0.013274787,    0.0146279344,    0.0161184594,    0.0177602228, 
	   0.0195684694,    0.0215599611,    0.0237531345,    0.0261682533, 
	   0.0288275983,     0.031755656,    0.0349793285,    0.0385281816, 
	   0.0424346849,    0.0467345081,    0.0514668077,    0.0566745736, 
	   0.0624049902,     0.068709828,    0.0756458789,    0.0832754225, 
	   0.0916667506,      0.10089469,      0.11104124,      0.12219622, 
	    0.134457946,     0.147934079,     0.162742332,     0.179011479, 
	    0.196882278,     0.216508493,     0.238058075,      0.26171428, 
	    0.287677079,     0.316164374,     0.347413689,     0.381683528, 
	    0.419255197,     0.460434496,     0.505553663,     0.554973185, 
	    0.609084129,     0.668310046,     0.733109355,     0.803977668, 
	    0.881450176,     0.966103852,      1.05856049,      1.15948844, 
	     1.26960564,      1.38968182,      1.52054048,      1.66306162, 
	     1.81818318,      1.98690259,      2.17027831,      2.36943054, 
	     2.58554149,       2.8198545,      3.07367301,      3.34835863, 
	     3.64532685,      3.96604347,      4.31201696,      4.68479109, 
	     5.08593607,      5.51703358,      5.97966528,       6.4753933, 
	     7.00574017,      7.57216644,      8.17604256,      8.81861877, 
	     9.50098896,      10.2240572,      10.9884892,      11.7946701, 
	     12.6426535,      13.5321083,      14.4622641,      15.4318485, 
	     16.4390354,      17.4813766,      18.5557499,      19.6583023, 
	     20.7844048,      21.9286003,      23.0845871,      24.2451916, 
	     25.4023743,      26.5472527,      27.6701469,      28.7606506, 
	     29.8077545,      30.7999763,      31.7255516,      32.5726585, 
	     33.3296852,      33.9855423,      34.5299873,      34.9540405, 
	     35.2503548,      35.4136658,      35.4411888,      35.3330536, 
	     35.0926628,      34.7270164,      34.2469559,      33.6672783, 
	     33.0067596,      32.2879753,      31.5369835,      30.7827816, 
	     30.0565758,      29.3908443,      28.8182049,      28.3701077, 
	     28.0754147,      27.9588795,      28.0396347,      28.3297291, 
	     28.8328266,      29.5431557,      30.4447823,      31.5113392, 
	     32.7062187,      33.9833488,      35.2885513,        36.56147, 
	     37.7380371,      38.7534256,      39.5453148,      40.0573654, 
	     40.2426643,      40.0669861,      39.5115738,       38.575264, 
	     37.2757378,      35.6496849,      33.7518463,      31.6527672, 
	     29.4354305,      27.1907902,      25.0125446,      22.9913769, 
	     21.2091446,      19.7334194,      18.6128502,      17.8737583, 
	     17.5182991,      17.5244007,      17.8475323,      18.4241753, 
	     19.1767006,      20.0192299,      20.8639565,      21.6273537, 
	     22.2357216,      22.6296597,      22.7671013,      22.6247749, 
	     22.1981106,      21.4997349,      20.5568256,      19.4077072, 
	     18.0980473,      16.6770229,      15.1937828,      13.6944523, 
	     12.2198143,      10.8037729,      9.47254848,      8.24456215, 
	     7.13086939,      6.13601923,      5.25915861,      4.49526453, 
	     3.83636904,      3.27267718,      2.79352903,      2.38815498, 
	     2.04623151,      1.75824165,      1.51567149,      1.31108141, 
	     1.13808215,     0.991257131,     0.866054714,     0.758673429, 
	    0.665954053,      0.58528477,     0.514521897,     0.451925486, 
	    0.396105468,      0.34597519,     0.300709218,     0.259702146, 
	    0.222527996,     0.188899443,     0.158628106,     0.131587014, 
	    0.107677303,    0.0868005008,    0.0688379183,    0.0536376387, 
	   0.0410092473,    0.0307254996,    0.0225295536,    0.0161460824, 
	   0.0112942234,   0.00770051079,   0.00511020282,   0.00329589355, 
	  0.00206282083,   0.00125087425,  0.000733686436,  0.000415528833, 
	 0.000226829725,  0.000119119861,    6.006084e-05,  2.90146272e-05, 
	 1.34002139e-05,  5.90308809e-06,  2.47442199e-06, 
/* Z=30 -> Zn */
	              0,   0.00758356694,   0.00835675746,   0.00920854881, 
	   0.0101468964,    0.0111805508,    0.0123191448,    0.0135732768, 
	   0.0149546033,    0.0164759494,    0.0181514192,    0.0199965239, 
	   0.0220283233,    0.0242655687,    0.0267288741,    0.0294409003, 
	    0.032426551,    0.0357131809,    0.0393308513,    0.0433125719, 
	   0.0476946011,    0.0525167435,    0.0578226894,    0.0636603907, 
	   0.0700824484,    0.0771465674,    0.0849160329,    0.0934602171, 
	    0.102855131,     0.113184072,     0.124538258,     0.137017578, 
	    0.150731325,     0.165799081,     0.182351604,     0.200531811, 
	    0.220495895,      0.24241434,     0.266473293,     0.292875737, 
	    0.321843058,     0.353616327,     0.388458163,      0.42665419, 
	    0.468515009,     0.514378071,     0.564609528,     0.619606674, 
	    0.679799855,     0.745654881,     0.817675412,     0.896405458, 
	     0.98243171,      1.07638621,      1.17894912,      1.29085076, 
	     1.41287434,       1.5458585,      1.69069958,      1.84835303, 
	     2.01983571,      2.20622754,      2.40867162,      2.62837481, 
	     2.86660838,      3.12470531,      3.40405941,      3.70612192, 
	     4.03239632,      4.38443279,      4.76382017,      5.17217684, 
	     5.61113739,      6.08233976,      6.58740616,      7.12792492, 
	     7.70542622,      8.32135487,      8.97703838,      9.67365646, 
	     10.4121981,      11.1934214,      12.0178061,      12.8855038, 
	     13.7962837,      14.7494745,       15.743906,      16.7778473, 
	     17.8489456,      18.9541626,      20.0897274,      21.2510662, 
	     22.4327774,       23.628582,      24.8313065,      26.0328789, 
	     27.2243538,      28.3959332,      29.5370655,      30.6365223, 
	     31.6825504,      32.6630554,       33.565815,      34.3787422, 
	     35.0902023,      35.6893539,      36.1665306,      36.5136604, 
	     36.7246971,      36.7960663,      36.7270737,      36.5203362, 
	     36.1821098,      35.7225685,      35.1559677,      34.5006943, 
	     33.7791481,      33.0174294,      32.2448692,      31.4933033, 
	     30.7961578,      30.1873283,      29.6998749,      29.3645439, 
	     29.2082367,      29.2524014,      29.5114975,      29.9915829, 
	     30.6891346,      31.5902023,      32.6699829,      33.8929062, 
	     35.2132759,      36.5765419,      37.9211655,      39.1810951, 
	     40.2887421,      41.1783295,      41.7895317,      42.0711174, 
	     41.9844475,      41.5065384,      40.6324844,      39.3769798, 
	     37.7747612,      35.8798294,      33.7633553,      31.5103321, 
	     29.2150574,      26.9756966,      24.8882465,      23.0403271, 
	     21.5052719,      20.3369961,      19.5661182,      19.1977043, 
	     19.2108955,      19.5604935,      20.1804256,      20.9887905, 
	     21.8940449,      22.8017826,      23.6214981,      24.2727242, 
	      24.690052,      24.8266621,      24.6561546,      24.1726494, 
	     23.3893356,      22.3357201,      21.0539761,      19.5947819, 
	     18.0130844,      16.3641319,      14.7000408,       13.067132, 
	     11.5040817,      10.0409126,      8.69875908,      7.49026728, 
	     6.42049551,      5.48813009,       4.6868639,      4.00679636, 
	     3.43573952,      2.96035767,      2.56709218,      2.24286032, 
	     1.97553825,      1.75425458,      1.56953263,      1.41331816, 
	     1.27893066,      1.16096592,      1.05517399,     0.958327234, 
	    0.868085802,     0.782866061,     0.701712489,     0.624174297, 
	    0.550186515,     0.479956865,     0.413860708,      0.35234639, 
	    0.295853883,     0.244749114,     0.199276358,     0.159528747, 
	    0.125437275,    0.0967761576,    0.0731818676,    0.0541823246, 
	   0.0392319411,    0.0277485177,    0.0191481132,    0.0128748761, 
	  0.00842383876,   0.00535572739,   0.00330392527,     0.001974585, 
	  0.00114143849,  0.000637121266,  0.000342774845,  0.000177418246, 
	  8.8172601e-05,   4.1987023e-05,  1.91159888e-05,  8.30204226e-06, 
	 3.43111788e-06,  1.34602396e-06,  4.99903308e-07, 
/* Z=31 -> Ga */
	              0,   0.00852765702,   0.00939538423,    0.0103511401, 
	   0.0114038112,    0.0125631774,    0.0138399927,    0.0152460942, 
	   0.0167944934,    0.0184995029,    0.0203768611,    0.0224438738, 
	   0.0247195587,    0.0272248294,    0.0299826581,    0.0330182984, 
	   0.0363594852,    0.0400366969,    0.0440833867,    0.0485363081, 
	   0.0534357913,    0.0588261187,    0.0647558719,    0.0712783411, 
	   0.0784519985,    0.0863409266,    0.0950153992,     0.104552411, 
	    0.115036309,     0.126559496,     0.139223069,     0.153137743, 
	    0.168424577,      0.18521595,     0.203656614,     0.223904639, 
	    0.246132687,     0.270529181,     0.297299653,     0.326668203, 
	     0.35887897,     0.394197792,     0.432913929,     0.475341916, 
	    0.521823466,     0.572729468,     0.628462255,     0.689457893, 
	    0.756188214,      0.82916379,     0.908935785,     0.996099114, 
	     1.09129465,      1.19521224,      1.30859292,      1.43223155, 
	     1.56697965,      1.71374738,      1.87350595,      2.04728913, 
	     2.23619533,      2.44138789,      2.66409683,      2.90561724, 
	     3.16730976,      3.45059752,      3.75696421,       4.0879488, 
	     4.44514084,      4.83017111,      5.24470425,      5.69042587, 
	     6.16902924,      6.68219566,      7.23157883,      7.81877804, 
	     8.44531345,       9.1125946,      9.82188702,      10.5742712, 
	     11.3706026,      12.2114649,      13.0971174,      14.0274382, 
	     15.0018721,      16.0193634,      17.0783005,      18.1764431, 
	     19.3108673,      20.4779015,      21.6730728,      22.8910522, 
	     24.1256256,      25.3696556,      26.6150837,       27.852932, 
	     29.0733395,      30.2656345,       31.418417,      32.5197067, 
	     33.5570984,      34.5179977,      35.3898659,      36.1605339, 
	     36.8185387,      37.3535347,      37.7566719,      38.0210876, 
	     38.1423149,      38.1187744,      37.9521637,      37.6478615, 
	     37.2152252,      36.6677933,      36.0233841,      35.3040123, 
	     34.5356483,       33.747757,      32.9726257,      32.2444992, 
	     31.5984364,      31.0690384,      30.6889477,      30.4872684, 
	     30.4879093,      30.7079563,      31.1561527,      31.8316021, 
	     32.7227669,      33.8068962,      35.0499611,      36.4071579, 
	     37.8240509,       39.238369,      40.5824127,      41.7860374, 
	     42.7800941,      43.5001793,       43.890461,      43.9074287, 
	     43.5232506,      42.7285004,      41.5340157,      39.9716415, 
	     38.0937195,      35.9712067,       33.690403,      31.3484421, 
	     29.0477142,      26.8895741,      24.9677486,      23.3619576, 
	     22.1322498,      21.3145809,      20.9180527,      20.9241257, 
	     21.2879238,      21.9416008,       22.799469,      23.7644787, 
	     24.7354374,      25.6143532,      26.3132248,      26.7597027, 
	     26.9012203,      26.7072792,      26.1698704,      25.3021564, 
	     24.1356525,      22.7163563,      21.1002407,      19.3485489, 
	     17.5233383,      15.6835661,      13.8819456,      12.1627388, 
	     10.5604572,      9.09944725,      7.79423141,      6.65043068, 
	     5.66610003,      4.83328152,      4.13962412,      3.56993246, 
	     3.10754323,      2.73547864,      2.43734455,      2.19797802, 
	     2.00386906,      1.84338927,      1.70686793,      1.58655405, 
	     1.47649825,       1.3723805,       1.2713052,      1.17157912, 
	     1.07248294,     0.974046946,     0.876836419,     0.781755924, 
	    0.689876378,     0.602290452,     0.519999146,     0.443830609, 
	    0.374389559,     0.312035769,     0.256886154,     0.208835483, 
	    0.167590007,     0.132708296,     0.103644118,    0.0797877386, 
	   0.0605023578,    0.0451543592,    0.0331367068,    0.0238855854, 
	   0.0168911163,    0.0117031559,   0.00793317705,    0.0052531953, 
	  0.00339253829,   0.00213308423,   0.00130346103,  0.000772663741, 
	 0.000443458732,  0.000245938398,   0.00013152929,    6.769006e-05, 
	 3.34491815e-05,  1.58351522e-05,  7.16499653e-06, 
/* Z=32 -> Ge */
	              0,   0.00956309866,    0.0105341999,    0.0116036041, 
	    0.012781214,    0.0140779186,    0.0155057013,     0.017077731, 
	   0.0188084971,    0.0207139235,    0.0228115208,    0.0251205303, 
	   0.0276621059,    0.0304594841,    0.0335382037,    0.0369263142, 
	   0.0406546257,    0.0447569638,    0.0492704771,    0.0542359389, 
	   0.0596980937,     0.065706037,    0.0723136291,    0.0795799345, 
	   0.0875697061,    0.0963539258,     0.106010377,     0.116624273, 
	    0.128288895,     0.141106412,     0.155188575,      0.17065765, 
	    0.187647358,     0.206303835,     0.226786688,     0.249270275, 
	    0.273944885,        0.301018,     0.330715984,     0.363285333, 
	    0.398994505,     0.438135654,     0.481026441,     0.528011978, 
	    0.579466999,     0.635797977,     0.697445393,      0.76488632, 
	    0.838636518,     0.919253469,      1.00733852,      1.10354006, 
	     1.20855594,      1.32313585,      1.44808495,      1.58426535, 
	     1.73259926,      1.89407134,       2.0697298,      2.26068926, 
	     2.46813154,      2.69330621,      2.93753028,      3.20218897, 
	     3.48873234,      3.79867363,      4.13358545,      4.49509239, 
	     4.88486719,       5.3046174,      5.75607777,      6.24099398, 
	     6.76110697,      7.31813383,      7.91374445,      8.54953384, 
	     9.22699547,      9.94748116,      10.7121696,      11.5220194, 
	     12.3777199,      13.2796488,      14.2278051,      15.2217598, 
	     16.2605915,      17.3428154,      18.4663296,      19.6283455, 
	     20.8253155,      22.0528927,      23.3058586,      24.5780907, 
	     25.8625298,      27.1511574,      28.4350147,      29.7042122, 
	      30.947998,      32.1548462,      33.3125725,      34.4085159, 
	     35.4297295,      36.3632393,      37.1963654,      37.9170303, 
	     38.5141754,      38.9781723,      39.3012772,      39.4780998, 
	     39.5060577,      39.3858528,      39.1218491,      38.7224464, 
	     38.2003059,      37.5724945,      36.8604622,      36.0898438, 
	     35.2900505,      34.4936523,      33.7355309,      33.0518074, 
	     32.4785271,      32.0502052,      31.7981815,      31.7489357, 
	     31.9223652,      32.3301697,       32.974411,      33.8463821, 
	     34.9258499,      36.1808281,      37.5679245,      39.0333176, 
	     40.5144424,      41.9423256,      43.2445335,      44.3486481, 
	     45.1861229,      45.6962967,      45.8303566,      45.5550232, 
	     44.8555984,      43.7381897,      42.2308273,      40.3832664, 
	     38.2653694,      35.9640236,      33.5786705,      31.2156677, 
	     28.9817562,      26.9771118,      25.2884769,      23.9829292, 
	     23.1028404,      22.6625137,      22.6468811,      23.0124207, 
	     23.6903152,      24.5915928,      25.6138363,      26.6488571, 
	     27.5906372,      28.3428307,      28.8251896,      28.9783802, 
	     28.7668781,      28.1798077,      27.2298222,      25.9503078, 
	     24.3912926,      22.6145477,      20.6883907,      18.6826401, 
	     16.6640873,       14.692812,      12.8194723,      11.0836477, 
	     9.51318455,      8.12442589,      6.92317629,      5.90617275, 
	      5.0628767,      4.37738991,      3.83033085,      3.40054822, 
	     3.06658769,      2.80786514,      2.60553789,      2.44308591, 
	     2.30663776,      2.18508101,      2.07000089,      1.95549023, 
	     1.83786941,      1.71534836,      1.58765924,      1.45568454, 
	     1.32110023,      1.18604922,      1.05285943,     0.923811674, 
	    0.800964117,     0.686030746,     0.580312252,     0.484671324, 
	    0.399544358,     0.324980497,     0.260698467,     0.206152767, 
	    0.160602391,     0.123176686,    0.0929344669,    0.0689144433, 
	   0.0501761921,    0.0358314328,    0.0250664186,    0.0171563402, 
	   0.0114727039,   0.00748487888,   0.00475682272,   0.00294012576, 
	  0.00176442915,   0.00102629955,  0.000577547122,  0.000313851546, 
	  0.00016437372,  8.27989352e-05,  4.00292702e-05,  1.85322497e-05, 
	 8.19730849e-06,     3.45587e-06,  1.38512553e-06, 
/* Z=33 -> As */
	              0,    0.0106970547,    0.0117810154,    0.0129744587, 
	   0.0142883882,    0.0157349017,    0.0173273031,    0.0190802179, 
	   0.0210097227,    0.0231334921,    0.0254709497,    0.0280434396, 
	   0.0308744106,    0.0339896269,     0.037417382,    0.0411887541, 
	   0.0453378633,    0.0499021634,    0.0549227707,    0.0604447946, 
	    0.066517733,    0.0731958821,    0.0805387571,      0.08861164, 
	   0.0974860564,     0.107240386,     0.117960483,     0.129740342, 
	     0.14268288,     0.156900704,     0.172516987,     0.189666405, 
	    0.208496153,     0.229167074,     0.251854748,     0.276750803, 
	    0.304064363,     0.334023327,     0.366876066,     0.402893037, 
	    0.442368507,     0.485622585,     0.533002973,     0.584887326, 
	    0.641685307,     0.703840911,      0.77183497,      0.84618777, 
	    0.927461326,      1.01626241,      1.11324513,      1.21911383, 
	      1.3346256,      1.46059346,      1.59788859,       1.7474432, 
	     1.91025293,       2.0873785,      2.27994871,      2.48916078, 
	     2.71628118,       2.9626472,       3.2296648,      3.51880932, 
	     3.83162022,      4.16970062,      4.53470898,      4.92835522, 
	     5.35238886,      5.80859041,      6.29875803,      6.82468891, 
	     7.38816309,      7.99091864,       8.6346283,      9.32086658, 
	     10.0510778,       10.826539,      11.6483126,      12.5172081, 
	     13.4337187,      14.3979797,      15.4096956,      16.4680901, 
	     17.5718288,      18.7189598,      19.9068489,       21.132103, 
	     22.3905201,      23.6770248,      24.9856186,       26.309351, 
	     27.6402931,      28.9695301,      30.2871914,      31.5824909, 
	      32.843811,       34.058815,      35.2146111,      36.2979393, 
	     37.2954292,      38.1938972,      38.9806633,      39.6439743, 
	     40.1733856,      40.5602608,       40.798214,      40.8836174, 
	     40.8160591,      40.5987968,      40.2391205,       39.748661, 
	     39.1435623,      38.4445152,      37.6766129,      36.8690186, 
	     36.0543747,      35.2680397,      34.5470161,      33.9286957, 
	     33.4493828,      33.1426659,      33.0376434,      33.1571846, 
	     33.5161629,      34.1199226,      34.9629745,      36.0280647, 
	     37.2857704,      38.6946411,      40.2020226,      41.7455711, 
	     43.2554817,      44.6573601,      45.8757401,      46.8379593, 
	     47.4783707,      47.7425194,       47.591114,      47.0034599, 
	     45.9800606,      44.5441399,      42.7418556,      40.6409988, 
	     38.3281746,      35.9044609,      33.4797516,      31.1660748, 
	     29.0702839,      27.2867222,       25.890379,      24.9311733, 
	     24.4299126,      24.3763313,      24.7295208,      25.4207458, 
	     26.3585052,      27.4353943,      28.5361843,      29.5463886, 
	     30.3605328,      30.8894272,      31.0658112,      30.8479671, 
	     30.2211285,      29.1966763,      27.8094196,      26.1133423, 
	     24.1763172,      22.0743446,      19.8858452,      17.6864281, 
	     15.5445194,      13.5180521,      11.6523371,      9.97909164, 
	     8.51652718,      7.27033091,      6.23531103,      5.39748621, 
	      4.7363801,      4.22732782,      3.84362626,      3.55840993, 
	     3.34617591,      3.18392992,      3.05195475,      2.93422794, 
	     2.81853795,      2.69635224,      2.56249619,      2.41469955, 
	     2.25306392,      2.07949781,      1.89715981,      1.70994186, 
	     1.52201581,      1.33746016,      1.15997362,     0.992677689, 
	    0.838000119,     0.697631419,     0.572540164,     0.463032812, 
	    0.368844569,     0.289247364,     0.223163977,     0.169279978, 
	    0.126146287,    0.0922683477,    0.0661795288,    0.0464977436, 
	   0.0319656767,    0.0214757733,    0.0140818423,    0.0089993896, 
	  0.00559721515,   0.00338271121,   0.00198329031,   0.00112615211, 
	 0.000618192425,  0.000327459391,  0.000167052596,  8.19083216e-05, 
	  3.8517559e-05,  1.73332683e-05,  7.44704357e-06,  3.04728474e-06, 
	 1.18457558e-06,  4.36291373e-07,  1.51825333e-07, 
/* Z=34 -> Se */
	              0,    0.0119372932,    0.0131442985,    0.0144729353, 
	   0.0159353968,    0.0175450891,    0.0193167478,    0.0212665647, 
	    0.023412345,    0.0257736463,    0.0283719618,    0.0312309004, 
	   0.0343763977,    0.0378369316,    0.0416437835,    0.0458312817, 
	   0.0504371151,    0.0555026494,     0.061073266,    0.0671987459, 
	   0.0739337057,    0.0813379958,    0.0894772559,    0.0984234065, 
	     0.10825526,     0.119059108,      0.13092947,     0.143969774, 
	    0.158293247,     0.174023628,     0.191296324,     0.210259274, 
	    0.231074095,     0.253917247,     0.278981388,     0.306476593, 
	    0.336631954,     0.369697094,     0.405943871,     0.445668072, 
	    0.489191413,     0.536863565,     0.589064121,     0.646205068, 
	    0.708732963,     0.777131438,     0.851923883,     0.933675826, 
	     1.02299809,      1.12054908,      1.22703815,       1.3432281, 
	     1.46993816,      1.60804641,      1.75849342,      1.92228353, 
	     2.10048819,      2.29424763,      2.50477171,      2.73334241, 
	     2.98131371,      3.25011039,      3.54122901,      3.85623431, 
	     4.19675541,      4.56448317,      4.96116066,      5.38857746, 
	     5.84855843,      6.34294891,      6.87360287,      7.44236088, 
	     8.05102921,      8.70135498,      9.39499664,      10.1334896, 
	      10.918211,      11.7503338,       12.630785,      13.5601883, 
	     14.5388126,       15.566514,      16.6426678,      17.7661018, 
	     18.9350376,      20.1470051,      21.3987942,      22.6863728, 
	     24.0048294,      25.3483276,      26.7100582,      28.0822124, 
	     29.4559669,      30.8215065,      32.1680565,      33.4839478, 
	     34.7567368,      35.9733429,      37.1202393,      38.1836853, 
	     39.1500206,      40.0060005,      40.7391586,      41.3382568, 
	     41.7937317,      42.0981712,      42.2468376,      42.2381363, 
	     42.0741005,      41.7608032,      41.3086815,      40.7327843, 
	     40.0528526,      39.2932434,      38.4826393,      37.6535606, 
	     36.8416176,      36.0845184,      35.4208412,      34.8885841, 
	     34.5235138,      34.3573875,      34.4161263,      34.7179756, 
	     35.2718391,      36.0758018,      37.1160431,      38.3662148, 
	      39.787384,        41.32864,      42.9284286,      44.5166397, 
	     46.0173683,      47.3523903,      48.4450989,      49.2248039, 
	     49.6311607,      49.6184044,      49.1591606,      48.2474632, 
	     46.9007454,      45.1604805,      43.0912971,      40.7785072, 
	     38.3239784,       35.840538,      33.4451866,      31.2515278, 
	     29.3619385,      27.8601265,      26.8046837,      26.2242584, 
	     26.1148643,      26.4396534,      27.1312809,      28.0967312, 
	     29.2242565,      30.3917904,      31.4761124,      32.3619385, 
	     32.9501152,      33.1642113,       32.954998,      32.3025703, 
	     31.2159939,      29.7307854,      27.9045391,      25.8112793, 
	     23.5351162,      21.1637859,      18.7826118,      16.4693298, 
	     14.2900476,      12.2965269,      10.5247993,      8.99505806, 
	     7.71264791,      6.66990852,      5.84862137,      5.22278023, 
	     4.76143408,      4.43138313,       4.1995635,      4.03500366, 
	     3.91029024,      3.80253935,      3.69388771,       3.5715704, 
	     3.42765212,      3.25850177,       3.0640893,      2.84719229, 
	     2.61257672,      2.36621928,      2.11460829,      1.86416125, 
	     1.62076914,      1.38947725,      1.17429018,     0.978091836, 
	    0.802656949,     0.648735344,      0.51618439,     0.404131949, 
	    0.311150461,     0.235428765,     0.174931139,     0.127536133, 
	   0.0911508054,    0.0637990758,    0.0436841026,    0.0292266831, 
	   0.0190826375,     0.012142837,   0.00751986308,   0.00452540815, 
	   0.0026422774,   0.00149435073,  0.000817191729,   0.00043132162, 
	 0.000219308145,  0.000107205211,  5.02776493e-05,  2.25724998e-05, 
	 9.67904725e-06,  3.95448069e-06,   1.5355223e-06,  5.65180528e-07, 
	 1.96645132e-07,  6.44895337e-08,  1.98743635e-08, 
/* Z=35 -> Br */
	              0,    0.0132921264,    0.0146331051,    0.0161088984, 
	   0.0177329853,    0.0195201822,    0.0214867741,    0.0236506481, 
	   0.0260314662,     0.028650824,    0.0315324441,    0.0347023793, 
	   0.0381892398,    0.0420244336,    0.0462424457,    0.0508811139, 
	   0.0559819713,    0.0615905747,    0.0677568987,    0.0745357424, 
	   0.0819871947,    0.0901771039,    0.0991776586,     0.109067902, 
	    0.119934432,     0.131872043,     0.144984514,     0.159385353, 
	    0.175198764,     0.192560494,      0.21161893,     0.232536137, 
	    0.255489141,      0.28067106,     0.308292657,     0.338583678, 
	    0.371794492,     0.408197671,     0.448090017,     0.491794199, 
	     0.53966099,     0.592071295,     0.649438441,     0.712210596, 
	    0.780873239,     0.855951786,      0.93801415,      1.02767372, 
	     1.12559199,      1.23248196,      1.34911025,       1.4763006, 
	     1.61493695,      1.76596534,      1.93039775,      2.10931349, 
	     2.30386162,      2.51526356,      2.74481344,      2.99387908, 
	     3.26390195,       3.5563972,      3.87295103,      4.21521711, 
	     4.58491421,       4.9838171,      5.41375208,      5.87658453, 
	     6.37420702,      6.90852499,      7.48143864,      8.09482002, 
	     8.75049114,      9.45019245,      10.1955519,      10.9880457, 
	     11.8289614,      12.7193441,      13.6599522,      14.6511993, 
	     15.6930943,      16.7851753,      17.9264507,      19.1153183, 
	     20.3495064,      21.6259899,      22.9409332,      24.2896233, 
	     25.6664085,      27.0646553,      28.4767075,      29.8938789, 
	     31.3064423,      32.7036819,      34.0739288,      35.4046707, 
	     36.6826897,      37.8942223,      39.0252113,      40.0615616, 
	     40.9894714,      41.7958069,      42.4685249,      42.9971275, 
	     43.3731651,      43.5907364,       43.647007,      43.5426941, 
	     43.2825279,      42.8756256,      42.3357697,      41.6815605, 
	     40.9364128,       40.128315,      39.2894249,      38.4553833, 
	     37.6643944,      36.9560318,      36.3698311,      35.9436417, 
	     35.7118492,      35.7034721,      35.9402695,      36.4349174, 
	     37.1894226,      38.1938171,      39.4253311,       40.848114, 
	       42.41362,      44.0617104,      45.7225609,      47.3192902, 
	     48.7713394,      49.9984169,      50.9248772,      51.4843483, 
	     51.6242218,       51.309845,       50.527977,      49.2892532, 
	     47.6293411,      45.6085396,      43.3097267,      40.8345222, 
	     38.2978897,       35.821312,      33.5249863,      31.5196171, 
	     29.8983116,      28.7293987,      28.0507507,      27.8662281, 
	     28.1446857,      28.8216991,      29.8040543,      30.9765835, 
	     32.2108383,      33.3748245,      34.3428879,      35.0049133, 
	     35.2739792,      35.0918694,      34.4320297,      33.2998734, 
	     31.7305603,      29.7846355,      27.5420494,      25.0952187, 
	     22.5417709,      19.9776058,      17.4907646,      15.1565361, 
	     13.0340099,      11.1641788,      9.56955624,      8.25513554, 
	      7.2104497,      6.41243458,      5.82876015,      5.42131615, 
	     5.14956951,      4.97355556,       4.8563447,      4.76587725, 
	     4.67614365,      4.56773186,      4.42781687,      4.24969149, 
	     4.03196096,      3.77752399,      3.49245358,      3.18488503, 
	     2.86398768,      2.53907728,      2.21890688,      1.91114998, 
	     1.62207043,       1.3563602,      1.11712253,     0.905964971, 
	     0.72317332,     0.567933977,     0.438579828,     0.332836568, 
	    0.248053551,     0.181407914,     0.130074203,     0.091356881, 
	   0.0627857521,    0.0421763659,    0.0276597235,    0.0176865794, 
	   0.0110120801,   0.00666654622,   0.00391813368,   0.00223209057, 
	  0.00123046106,  0.000655213196,  0.000336397381,  0.000166201702, 
	 7.88584803e-05,  3.58563702e-05,  1.55890593e-05,  6.46537183e-06, 
	 2.55166083e-06,  9.55858468e-07,   3.3895185e-07,  1.13457773e-07, 
	 3.57440975e-08,  1.05660556e-08,  2.92123437e-09, 
/* Z=36 -> Kr */
	              0,    0.0147705059,    0.0162571836,    0.0178929605, 
	   0.0196927078,    0.0216727629,    0.0238510743,    0.0262473654, 
	   0.0288832951,    0.0317826569,    0.0349715725,    0.0384787247, 
	   0.0423356034,    0.0465767719,    0.0512401536,    0.0563673638, 
	   0.0620040484,    0.0682002753,    0.0750109479,    0.0824962482, 
	   0.0907221287,    0.0997608677,     0.109691642,     0.120601147, 
	    0.132584304,     0.145744994,     0.160196885,     0.176064253, 
	    0.193483025,     0.212601721,     0.233582556,     0.256602705, 
	    0.281855494,     0.309551865,     0.339921743,     0.373215765, 
	    0.409706801,     0.449691921,     0.493494183,     0.541464806, 
	     0.59398514,     0.651469171,     0.714365721,     0.783161163, 
	    0.858381927,     0.940597057,       1.0304215,      1.12851858, 
	     1.23560321,      1.35244453,      1.47986948,      1.61876547, 
	     1.77008295,      1.93483865,       2.1141181,      2.30907798, 
	     2.52094722,      2.75103045,      3.00070715,       3.2714324, 
	     3.56473708,      3.88222599,      4.22557497,      4.59652662, 
	     4.99688578,      5.42851114,      5.89330673,       6.3932085, 
	     6.93017197,      7.50615311,      8.12308884,      8.78287315, 
	     9.48732758,      10.2381706,      11.0369816,      11.8851566, 
	     12.7838678,      13.7340097,      14.7361422,      15.7904348, 
	     16.8966026,      18.0538311,      19.2607212,      20.5151997, 
	     21.8144569,      23.1548748,      24.5319576,      25.9402637, 
	     27.3733692,      28.8238068,      30.2830582,      31.7415371, 
	     33.1886215,      34.6126862,      36.0012054,      37.3408623, 
	     38.6177216,      39.8174477,      40.9255447,         41.9277, 
	     42.8101387,      43.5600281,      44.1659584,       44.618412, 
	      44.910305,      45.0375061,      44.9993439,      44.7991104, 
	     44.4444695,      43.9477844,      43.3263397,      42.6023865, 
	     41.8030052,      40.9597511,      40.1080513,      39.2863617, 
	     38.5350342,      37.8949509,      37.4059029,      37.1048088, 
	     37.0237885,      37.1881752,      37.6146202,      38.3093147, 
	     39.2665291,      40.4675598,      41.8802147,      43.4589272, 
	     45.1456528,      46.8715096,      48.5592613,      50.1265335, 
	     51.4897461,      52.5685196,      53.2904396,      53.5958176, 
	     53.4422264,      52.8083992,      51.6972198,      50.1374054, 
	     48.1836891,      45.9152451,      43.4322739,      40.8508568, 
	     38.2962341,      35.8948631,      33.7658081,      32.0121231, 
	     30.7128391,      29.9164047,      29.6361732,      29.8484726, 
	     30.4935913,      31.4796867,      32.6894112,      33.9887085, 
	     35.2370071,      36.2979202,      37.0494576,      37.3928757, 
	     37.2593117,      36.6138382,      35.4565506,      33.8208466, 
	     31.7691936,      29.3869152,      26.7746868,      24.0404415, 
	     21.2914333,      18.6270638,      16.1329403,      13.8765774, 
	     11.9048424,      10.2431936,       8.8965559,      7.85158587, 
	     7.07998753,      6.54248142,      6.19303942,       5.9830122, 
	     5.86483717,      5.79509354,      5.73676443,      5.66064692, 
	     5.54593754,      5.38009739,      5.15812206,      4.88139343, 
	     4.55628014,      4.19264507,       3.8023994,      3.39820099, 
	     2.99237514,      2.59608626,       2.2187736,      1.86783135, 
	     1.54850328,      1.26394939,      1.01544094,     0.802638113, 
	    0.623913646,     0.476687282,     0.357747257,     0.263539672, 
	    0.190414235,     0.134820879,    0.0934552029,    0.0633558258, 
	   0.0419584028,    0.0271130092,    0.0170728061,    0.0104618277, 
	  0.00622956781,   0.00359909236,   0.00201425515,   0.00109015417, 
	 0.000569563126,  0.000286725786,  0.000138807649,  6.44898537e-05, 
	 2.86924333e-05,  1.21972798e-05,  4.94256983e-06,  1.90441835e-06, 
	 6.95930453e-07,  2.40538498e-07,  7.84124552e-08,  2.40368365e-08, 
	 6.90739466e-09,  1.85477467e-09,  4.63815708e-10, 
/* Z=37 -> Rb */
	              0,     0.016384406,    0.0180295967,    0.0198393706, 
	   0.0218301006,    0.0240197647,    0.0264281165,    0.0290768482, 
	   0.0319897793,    0.0351930633,    0.0387154147,    0.0425883457, 
	   0.0468464456,    0.0515276603,    0.0566736236,    0.0623300038, 
	   0.0685468763,    0.0753791481,    0.0828870013,    0.0911363885, 
	     0.10019958,     0.110155709,      0.12109147,     0.133101717, 
	    0.146290302,     0.160770789,     0.176667437,     0.194116026, 
	    0.213264972,     0.234276414,     0.257327348,     0.282610983, 
	     0.31033805,     0.340738356,      0.37406233,     0.410582662, 
	    0.450596213,     0.494425863,     0.542422593,     0.594967723, 
	    0.652475059,     0.715393364,     0.784209073,     0.859448552, 
	    0.941681385,      1.03152275,      1.12963641,      1.23673844, 
	     1.35359907,      1.48104692,      1.61997116,      1.77132511, 
	     1.93612885,      2.11547232,       2.3105166,      2.52249813, 
	     2.75272822,       3.0025959,      3.27356791,      3.56718802, 
	     3.88507724,         4.22893,      4.60051107,      5.00165176, 
	     5.43424082,      5.90021706,      6.40155792,      6.94026518, 
	     7.51834965,      8.13780975,      8.80060863,      9.50864792, 
	     10.2637358,      11.0675516,      11.9216061,      12.8271923, 
	     13.7853441,       14.796773,      15.8618135,        16.98036, 
	     18.1517925,      19.3749123,      20.6478653,       21.968071, 
	     23.3321419,      24.7358208,      26.1739082,      27.6402073, 
	     29.1274757,      30.6273937,      32.1305504,       33.626461, 
	     35.1035919,      36.5494385,      37.9506454,      39.2931404, 
	     40.5623512,      41.7434425,      42.8216209,      43.7824974, 
	     44.6124878,      45.2992668,      45.8322601,      46.2031822, 
	      46.406559,      46.4402771,      46.3061028,      46.0101204, 
	     45.5631523,      44.9809837,      44.2845192,      43.4997063, 
	     42.6572647,      41.7921677,      40.9428825,      40.1502991, 
	     39.4564362,       38.902874,      38.5289497,      38.3698311, 
	     38.4544907,       38.803688,      39.4280548,      40.3264656, 
	     41.4847488,      42.8749237,      44.4550667,      46.1699028, 
	     47.9522095,      49.7250671,      51.4048958,      52.9053001, 
	     54.1414795,      55.0350647,      55.5191307,       55.543045, 
	     55.0768166,      54.1146126,      52.6770477,      50.8119621, 
	     48.5934296,      46.1188889,      43.5043335,      40.8778114, 
	     38.3714333,      36.1125069,      34.2143517,      32.7675476, 
	     31.8324318,      31.4335213,      31.5565434,      32.1484718, 
	     33.1207199,      34.3553429,      35.7137947,      37.0474815, 
	     38.2092171,      39.0644684,      39.5014496,      39.4390793, 
	     38.8322487,      37.6739235,      35.9941063,      33.8558884, 
	     31.3490906,      28.5822468,      25.6736469,      22.7423134, 
	     19.8996277,      17.2422371,       14.846715,      12.7662144, 
	     11.0292349,      9.64038086,      8.58285141,      7.82228708, 
	     7.31151199,      6.99567366,       6.8173089,       6.7209177, 
	     6.65673113,      6.58345604,      6.46992779,      6.29568291, 
	     6.05059004,      5.73372078,      5.35168791,      4.91668797, 
	     4.44446611,      3.95238495,      3.45773935,      2.97640157, 
	     2.52184391,      2.10453963,      1.73171222,       1.4073801, 
	       1.132635,     0.906080246,     0.724364817,      0.58274883, 
	    0.475652605,     0.397148907,     0.341374457,     0.302846611, 
	    0.276684344,     0.258738846,     0.245647579,     0.234826341, 
	    0.224417076,     0.213205993,     0.200525284,     0.186148554, 
	    0.170186684,     0.152989283,     0.135054633,     0.116950028, 
	   0.0992445722,    0.0824552476,    0.0670074001,    0.0532102846, 
	   0.0412470363,    0.0311778542,    0.0229541399,    0.0164404195, 
	   0.0114405444,   0.00772457942,    0.0050533074,   0.00319810049, 
	  0.00195492571,   0.00115227827,  0.000653731928, 
/* Z=38 -> Sr */
	              0,    0.0181436762,    0.0199610479,     0.021959763, 
	   0.0241578054,    0.0265749283,    0.0292328242,    0.0321553051, 
	    0.035368517,     0.038901154,    0.0427847095,    0.0470537432, 
	    0.051746171,    0.0569035783,    0.0625715703,    0.0688001737, 
	   0.0756442025,    0.0831637457,    0.0914246589,     0.100499056, 
	    0.110465974,      0.12141189,     0.133431494,     0.146628395, 
	    0.161115929,     0.177018076,     0.194470316,     0.213620707, 
	    0.234630987,     0.257677764,     0.282953739,     0.310669214, 
	    0.341053367,     0.374356061,     0.410849363,     0.450829417, 
	    0.494618386,     0.542566597,     0.595054448,     0.652494967, 
	    0.715336204,     0.784063697,       0.8592031,     0.941323221, 
	     1.03103864,      1.12901282,      1.23596096,      1.35265374, 
	     1.47991979,      1.61864865,      1.76979494,      1.93438041, 
	     2.11349678,      2.30830884,      2.52005696,      2.75005817, 
	     2.99970841,      3.27048302,      3.56393671,      3.88170314, 
	     4.22549343,       4.5970912,       4.9983511,      5.43118811, 
	     5.89757442,      6.39952421,      6.93908358,      7.51831436, 
	     8.13927269,      8.80399132,      9.51444626,      10.2725353, 
	     11.0800352,      11.9385681,      12.8495541,      13.8141623, 
	     14.8332596,      15.9073486,      17.0365067,      18.2203102, 
	     19.4577808,      20.7472858,      22.0864868,      23.4722481, 
	     24.9005661,      26.3665028,      27.8641205,      29.3864288, 
	     30.9253407,      32.4716606,      34.0150719,      35.5441856, 
	      37.046566,      38.5088654,      39.9169273,      41.2560005, 
	     42.5109482,      43.6665573,      44.7078705,       45.620575, 
	      46.391468,      47.0089302,      47.4634666,      47.7482491, 
	     47.8596687,      47.7978935,      47.5673447,      47.1771545, 
	     46.6414642,       45.979641,      45.2162704,      44.3810005, 
	     43.5080948,      42.6357422,      41.8051262,      41.0591469, 
	     40.4409294,      39.9920692,       39.750679,      39.7493362, 
	       40.01297,      40.5568466,      41.3847351,      42.4874191, 
	     43.8416824,      45.4098587,      47.1401787,      48.9678307, 
	     50.8169746,      52.6035614,       54.239006,      55.6345482, 
	     56.7061462,      57.3796616,      57.5959892,      57.3158646, 
	     56.5238876,      55.2314377,      53.4781303,      51.3314819, 
	     48.8846474,      46.2521248,      43.5635414,      40.9557877, 
	     38.5639534,      36.5117035,      34.9017944,      33.8076134, 
	     33.2664871,       33.275547,      33.7906456,      34.7286263, 
	     35.9729004,      37.3819466,       38.800106,      40.0696182, 
	      41.042942,      41.5940895,      41.6280708,      41.0875206, 
	     39.9560814,      38.2582359,      36.0558662,      33.4419479, 
	     30.5320683,      27.4546623,      24.3408566,      21.3147945, 
	       18.48522,      15.9389238,       13.736454,      11.9102602, 
	     10.4652357,      9.38140297,      8.61832047,      8.12067604, 
	     7.82445002,      7.66306543,      7.57295513,      7.49814272, 
	     7.39353895,       7.2268281,       6.9789834,      6.64356279, 
	     6.22503757,      5.73646736,      5.19682932,      4.62830639, 
	     4.05377293,       3.4946723,      2.96939373,      2.49220228, 
	      2.0727129,      1.71585464,       1.4222368,      1.18881059, 
	     1.00971508,     0.877196789,     0.782509506,     0.716719747, 
	    0.671366751,     0.638949931,     0.613238871,     0.589419246, 
	     0.56409812,     0.535202682,     0.501803458,     0.463895261, 
	    0.422161162,     0.377741486,     0.332022369,     0.286454201, 
	    0.242405787,      0.20105648,     0.163326055,     0.129839987, 
	    0.100925893,    0.0766359642,    0.0567885265,    0.0410217904, 
	   0.0288528223,    0.0197351631,    0.0131097864,   0.00844580494, 
	  0.00526900776,   0.00317814341,   0.00185033865,   0.00103800255, 
	 0.000560028711,  0.000290027499,  0.000143877318, 
/* Z=39 -> Y */
	              0,    0.0200563129,    0.0220601913,    0.0242634956, 
	   0.0266859476,    0.0293492004,    0.0322770253,    0.0354955159, 
	   0.0390333086,    0.0429218262,    0.0471955426,    0.0518922769, 
	   0.0570535026,    0.0627247021,    0.0689557493,    0.0758012906, 
	   0.0833212286,    0.0915811956,      0.10065309,     0.110615619, 
	    0.121554993,      0.13356553,     0.146750465,     0.161222726, 
	    0.177105755,      0.19453451,     0.213656485,     0.234632745, 
	    0.257639199,     0.282867759,     0.310527802,     0.340847671, 
	    0.374076128,     0.410484135,     0.450366735,     0.494044781, 
	    0.541867197,     0.594213009,     0.651493669,     0.714155734, 
	    0.782682955,     0.857599497,     0.939472318,      1.02891433, 
	     1.12658715,      1.23320484,      1.34953618,       1.4764086, 
	     1.61471081,      1.76539648,      1.92948723,      2.10807514, 
	      2.3023262,      2.51348281,      2.74286461,      2.99187231, 
	     3.26198721,      3.55477238,      3.87187099,      4.21500731, 
	     4.58598089,      4.98666573,       5.4190011,      5.88498688, 
	     6.38667202,      6.92614269,      7.50550747,      8.12687778, 
	     8.79235077,      9.50397682,      10.2637358,      11.0735054, 
	     11.9350119,       12.849802,      13.8191833,      14.8441753, 
	     15.9254503,      17.0632725,      18.2574253,      19.5071449, 
	     20.8110428,      22.1670265,      23.5722275,      25.0229168, 
	     26.5144386,      28.0411377,      29.5963039,      31.1721191, 
	     32.7596359,      34.3487663,      35.9282875,      37.4859009, 
	     39.0083046,      40.4813347,      41.8901024,      43.2192421, 
	     44.4531708,      45.5764122,      46.5739899,      47.4318504, 
	     48.1373596,      48.6798248,      49.0510559,      49.2459259, 
	     49.2629395,      49.1047745,      48.7787476,      48.2971992, 
	     47.6777611,      46.9434471,      46.1225662,      45.2483711, 
	     44.3585052,      43.4940948,      42.6986275,      42.0164986, 
	     41.4913063,      41.1639252,      41.0704346,      41.2399406, 
	     41.6924591,      42.4369125,      43.4694366,      44.7721062, 
	     46.3122177,      48.0422783,      49.9007988,       51.813942, 
	     53.6981163,      55.4634247,      57.0178833,      58.2723389, 
	     59.1457138,      59.5704193,      59.4975357,      58.9013443, 
	     57.7828827,      56.1720772,      54.1281853,      51.7381973, 
	     49.1131897,      46.3825645,      43.6864471,      41.1665878, 
	       38.95644,      37.1710358,      35.8976212,       35.187851, 
	     35.0523872,      35.4585457,      36.3314323,      37.5586319, 
	     38.9982033,      40.4893951,      41.8651237,      42.9651375, 
	     43.6486282,      43.8051262,      43.3627129,      42.2928085, 
	     40.6112556,      38.3756409,      35.6793175,      32.6427689, 
	     29.4031906,      26.1033344,      22.8805408,      19.8569527, 
	     17.1316013,      14.7749882,      12.8263788,      11.2938967, 
	      10.157176,      9.37213421,      8.87726498,      8.60072613, 
	     8.46750355,      8.40597534,      8.35332108,      8.25942039, 
	      8.0890522,      7.82242727,      7.45425272,      6.99163246, 
	     6.45119429,      5.85585213,      5.23155689,      4.60436583, 
	     3.99802661,      3.43222523,      2.92152834,      2.47499084, 
	     2.09633827,       1.7845943,      1.53500891,      1.34013617, 
	     1.19093335,      1.07777035,     0.991275728,      0.92297554, 
	    0.865714252,     0.813868105,     0.763384223,     0.711682856, 
	    0.657467425,     0.600480437,     0.541240156,     0.480782866, 
	     0.42042914,     0.361584693,     0.305582374,     0.253565818, 
	    0.206414267,     0.164704666,     0.128706619,    0.0984039232, 
	   0.0735362172,    0.0536536872,    0.0381776541,    0.0264606923, 
	   0.0178408138,    0.0116857383,   0.00742503908,    0.0045696306, 
	  0.00271960627,   0.00156256941,  0.000865183189,  0.000460785843, 
	 0.000235589978,  0.000115393501,  5.40289911e-05, 
/* Z=40 -> Zr */
	              0,    0.0221357103,    0.0243416056,    0.0267664306, 
	   0.0294317771,    0.0323613398,    0.0355811268,    0.0391196758, 
	   0.0430082865,    0.0472813137,      0.05197642,    0.0571349151, 
	   0.0628020912,    0.0690275952,    0.0758658424,    0.0833764598, 
	    0.091624774,     0.100682311,      0.11062742,     0.121545874, 
	    0.133531496,     0.146687001,     0.161124706,     0.176967382, 
	    0.194349304,     0.213417113,     0.234330997,     0.257265896, 
	    0.282412648,     0.309979469,     0.340193421,     0.373301893, 
	     0.40957433,     0.449304104,     0.492810369,     0.540440142, 
	    0.592570484,     0.649610698,     0.712005019,     0.780234873, 
	    0.854821861,     0.936330378,      1.02537072,      1.12260175, 
	     1.22873485,      1.34453607,       1.4708302,      1.60850358, 
	     1.75850761,      1.92186201,      2.09965754,      2.29305935, 
	     2.50330949,      2.73172879,      2.97972012,      3.24876785, 
	     3.54044032,      3.85638905,      4.19834709,      4.56812716, 
	     4.96761942,      5.39878368,      5.86364508,      6.36428261, 
	     6.90281868,      7.48140574,      8.10220718,      8.76737881, 
	     9.47904301,      10.2392607,      11.0500011,      11.9131012, 
	      12.830225,      13.8028183,      14.8320532,      15.9187727, 
	     17.0634308,       18.266016,      19.5259914,      20.8422127, 
	     22.2128468,      23.6353016,      25.1061382,      26.6210041, 
	     28.1745472,      29.7603607,      31.3709316,      32.9975929, 
	     34.6305122,      36.2587013,      37.8700256,      39.4513054, 
	     40.9883957,      42.4663582,      43.8696442,      45.1823692, 
	     46.3885994,      47.4727402,      48.4199486,      49.2166176, 
	     49.8508873,      50.3132095,      50.5969276,      50.6988564, 
	     50.6198578,      50.3653412,      49.9457169,       49.376709, 
	      48.679554,        47.88097,      47.0129433,      46.1122398, 
	     45.2196312,      44.3788643,      43.6352806,      43.0341949, 
	     42.6190071,      42.4291229,      42.4977493,      42.8496704, 
	     43.4991264,      44.4479065,      45.6838264,      47.1797409, 
	     48.8931618,      50.7667542,      52.7296219,      54.6995659, 
	     56.5862999,      58.2954674,      59.7334595,      60.8127327, 
	     61.4573936,      61.6086922,       61.230072,      60.3113213, 
	     58.8714447,      56.9598579,      54.6556244,      52.0645065, 
	     49.3138428,      46.5453186,      43.9060478,      41.5384445, 
	     39.5696106,      38.1011505,      37.2002258,      36.8928642, 
	     37.1602058,      37.9383125,      39.1217155,       40.570652, 
	     42.1214294,      43.5990982,       44.831337,      45.6621971, 
	     45.9645157,      45.6497955,      44.6746635,       43.043396, 
	     40.8063622,      38.0546722,      34.9116058,      31.5217896, 
	     28.0391026,      24.6144695,      21.3845844,      18.4625473, 
	     15.9310389,      13.8385553,      12.1987963,      10.9930611, 
	     10.1752167,      9.67856216,      9.42379665,      9.32722664, 
	     9.30837727,      9.29636002,      9.23449898,      9.08298016, 
	     8.81954384,      8.43842697,      7.94793081,      7.36708641, 
	     6.72189283,      6.04158354,      5.35527945,      4.68929434, 
	     4.06522131,      3.49885154,      2.99984527,      2.57205224, 
	     2.21430516,      1.92150831,      1.68585634,      1.49802744, 
	     1.34824383,      1.22712755,      1.12631583,      1.03883827, 
	    0.959281027,     0.883778036,     0.809877396,     0.736328065, 
	    0.662828028,     0.589764833,     0.517971456,      0.44851172, 
	    0.382502854,     0.320977837,     0.264787287,     0.214538515, 
	    0.170567349,     0.132938504,     0.101468757,    0.0757670552, 
	   0.0552853309,     0.039373979,     0.027336264,    0.0184769854, 
	   0.0121419011,   0.00774597656,   0.00479002018,   0.00286664558, 
	  0.00165750389,  0.000924294232,  0.000496173277,  0.000255902996, 
	 0.000126544939,  5.98694314e-05,  2.70378823e-05, 
/* Z=41 -> Nb */
	              0,    0.0243913997,    0.0268155895,    0.0294797029, 
	   0.0324073322,    0.0356243588,     0.039159175,    0.0430429392, 
	   0.0473098159,    0.0519972704,    0.0571463853,     0.062802203, 
	   0.0690140724,    0.0758361071,    0.0833275542,    0.0915533379, 
	    0.100584559,     0.110499062,     0.121382028,     0.133326724, 
	    0.146435127,     0.160818815,     0.176599726,     0.193911195, 
	    0.212898865,     0.233721778,      0.25655365,     0.281583935, 
	    0.309019446,     0.339085549,     0.372027874,     0.408114046, 
	    0.447635293,     0.490908563,     0.538278401,     0.590119302, 
	     0.64683789,     0.708875299,     0.776709974,     0.850860059, 
	    0.931886554,      1.02039599,      1.11704385,      1.22253704, 
	     1.33763826,      1.46316803,      1.60000873,       1.7491082, 
	     1.91148245,      2.08821893,      2.28048015,      2.48950553, 
	      2.7166152,      2.96321106,      3.23077869,      3.52088904, 
	     3.83519721,      4.17544413,      4.54345131,      4.94112158, 
	     5.37043238,      5.83342838,      6.33221626,      6.86895084, 
	     7.44582319,       8.0650444,      8.72882652,       9.4393568, 
	     10.1987753,      11.0091391,      11.8723888,       12.790307, 
	      13.764473,      14.7962103,      15.8865318,      17.0360794, 
	     18.2450523,      19.5131397,       20.839447,      22.2224121, 
	     23.6597271,      25.1482601,      26.6839695,      28.2618313, 
	     29.8757725,      31.5186043,      33.1819839,      34.8563728, 
	     36.5310516,      38.1941109,      39.8325424,      41.4322853, 
	     42.9784088,      44.4552574,      45.8466988,      47.1364288, 
	     48.3082962,      49.3467293,      50.2372017,      50.9667244, 
	     51.5244293,      51.9021339,      52.0949478,      52.1018639, 
	     51.9263153,      51.5766525,      51.0665436,      50.4152336, 
	     49.6476212,      48.7941284,      47.8903313,      46.9763031, 
	     46.0956726,      45.2943726,      44.6190834,      44.1154099, 
	     43.8258438,      43.7875137,      44.0299301,      44.5727234, 
	     45.4235611,       46.576416,      48.0102386,      49.6883087, 
	     51.5583229,      53.5533257,      55.5936127,      57.5896034, 
	     59.4456291,      61.0645866,      62.3531914,      63.2277069, 
	     63.6196747,      63.4814072,      62.7907181,      61.5545006, 
	     59.8107605,       57.628685,      55.1065521,      52.3673553, 
	     49.5521431,      46.8114433,      44.2951355,      42.1415405, 
	     40.4664955,      39.3534164,      38.8452835,      38.9394188, 
	     39.5857964,      40.6892166,      42.1154633,      43.7010078, 
	     45.2655525,      46.6263123,      47.6127663,       48.080452, 
	     47.9225655,       47.078228,      45.5367393,      43.3374748, 
	     40.5655327,       37.343647,      33.8212738,      30.1618347, 
	     26.5294437,      23.0762138,      19.9313087,      17.1925831, 
	     14.9214563,      13.1412373,      11.8388929,      10.9698334, 
	      10.465024,      10.2395754,      10.2018137,      10.2619104, 
	     10.3392296,      10.3678436,      10.2998505,      10.1064987, 
	     9.77732658,      9.31770039,      8.74531746,      8.08617783, 
	     7.37056208,      6.62941313,      5.89140511,      5.18086052, 
	     4.51652861,      3.91117048,       3.3717916,      2.90035939, 
	     2.49480772,      2.15016055,      1.85962808,      1.61557746, 
	     1.41030991,      1.23662591,      1.08818007,     0.959656298, 
	    0.846800327,     0.746353209,     0.655924261,     0.573837101, 
	    0.498970896,     0.430613995,     0.368336737,     0.311887085, 
	    0.261107981,     0.215875104,     0.176051736,     0.141458571, 
	    0.111856401,    0.0869397745,    0.0663399473,    0.0496352874, 
	   0.0363672227,    0.0260592792,    0.0182371754,    0.0124476021, 
	  0.00827398151,   0.00534795178,   0.00335600157,   0.00204130239, 
	  0.00120143581,  0.000683011778,   0.00037435125,  0.000197427435, 
	 9.99829499e-05,  4.85182754e-05,  2.25098047e-05, 
/* Z=42 -> Mo */
	              0,    0.0268441532,    0.0295047909,    0.0324280001, 
	    0.035639517,    0.0391675718,    0.0430431254,    0.0473001339, 
	   0.0519758314,    0.0571110286,    0.0627504662,    0.0689431652, 
	   0.0757428557,    0.0832083747,    0.0914041698,     0.100400843, 
	    0.110275619,     0.121113077,     0.133005753,     0.146054834, 
	     0.16037102,      0.17607531,     0.193299934,     0.212189376, 
	    0.232901394,     0.255608261,     0.280497909,     0.307775438, 
	    0.337664306,     0.370408237,     0.406272501,     0.445546031, 
	    0.488543153,     0.535605669,     0.587104976,     0.643444479, 
	    0.705061853,      0.77243185,      0.84606874,     0.926529348, 
	     1.01441586,      1.11037922,      1.21512163,      1.32940066, 
	     1.45403218,      1.58989358,      1.73792756,      1.89914501, 
	     2.07462907,       2.2655375,      2.47310638,      2.69865227, 
	     2.94357467,      3.20935822,      3.49757385,      3.80987906, 
	     4.14801788,      4.51382017,      4.90919733,      5.33614063, 
	     5.79671288,      6.29304361,      6.82731581,      7.40175724, 
	     8.01862144,      8.68017197,      9.38866043,      10.1462984, 
	     10.9552279,      11.8174896,      12.7349806,      13.7094107, 
	     14.7422523,      15.8346863,      16.9875393,      18.2012234, 
	     19.4756546,      20.8101902,      22.2035427,      23.6536942, 
	     25.1578274,      26.7122307,      28.3122215,      29.9520721, 
	     31.6249542,      33.3228645,      35.0365982,      36.7557411, 
	     38.4686584,      40.1625519,       41.823513,      43.4366646, 
	     44.9863091,      46.4561501,      47.8295593,      49.0899162, 
	     50.2209892,      51.2073784,      52.0350456,      52.6918259, 
	     53.1680489,      53.4571266,      53.5561714,      53.4665756, 
	     53.1945686,      52.7516365,      52.1548805,       51.427166, 
	     50.5971146,      49.6988258,      48.7713547,      47.8578796, 
	     47.0045509,      46.2590675,      45.6689148,      45.2793465, 
	     45.1311989,      45.2585602,      45.6864052,      46.4283867, 
	     47.4848518,      48.8412704,      50.4672318,      52.3161659, 
	     54.3258896,      56.4201126,      58.5109482,      60.5023804, 
	      62.294693,      63.7896385,      64.8961639,      65.5363922, 
	     65.6514435,      65.2067337,      64.1962433,      62.6454048, 
	     60.6120605,      58.1853104,      55.4819679,      52.6406593, 
	     49.8136978,      47.1571579,      44.8197708,      42.9314346, 
	     41.5922966,      40.8634529,      40.7602081,      41.2487259, 
	     42.2466621,      43.6279793,        45.23172,       46.874157, 
	     48.3633041,      49.5144615,      50.1653786,      50.1896172, 
	     49.5067787,      48.0887184,      45.9611626,      43.2006645, 
	      39.927269,      36.2936821,      32.4720459,      28.6395588, 
	     24.9643135,      21.5926228,      18.6388798,      16.1788254, 
	     14.2465954,      12.8356562,      11.9032927,      11.3779516, 
	     11.1685085,      11.1743708,      11.2953072,      11.4400291, 
	     11.5327883,      11.5175428,      11.3596096,      11.0449924, 
	      10.577836,      9.97660065,      9.26957417,      8.49032879, 
	     7.67359781,      6.85192776,      6.05326509,      5.29952955, 
	     4.60608625,      3.98196077,      3.43058991,      2.95089936, 
	     2.53851032,       2.1869216,      1.88855135,      1.63557506, 
	     1.42053366,      1.23672283,      1.07839191,     0.940796435, 
	    0.820149124,     0.713509083,     0.618643999,     0.533887982, 
	    0.458010107,     0.390100867,     0.329478741,     0.275615603, 
	    0.228078485,     0.186485246,     0.150471449,     0.119667351, 
	   0.0936834067,    0.0721039176,    0.0544875115,    0.0403736494, 
	   0.0292933956,    0.0207828488,    0.0143972905,   0.00972423237, 
	  0.00639391271,   0.00408630027,   0.00253417785,   0.00152247713, 
	 0.000884514069,  0.000496023451,  0.000267984811,  0.000139206924, 
	 6.93818947e-05,  3.31067822e-05,  1.50896549e-05, 
/* Z=43 -> Tc */
	              0,    0.0295106377,    0.0324273147,    0.0356309786, 
	   0.0391496755,    0.0430141538,    0.0472581275,    0.0519185588, 
	   0.0570359565,    0.0626547262,    0.0688235164,    0.0755956173, 
	   0.0830294192,    0.0911888406,      0.10014388,     0.109971151, 
	    0.120754488,      0.13258563,     0.145564899,     0.159802005, 
	    0.175416887,     0.192540556,     0.211316213,       0.2319002, 
	    0.254463166,     0.279191345,     0.306287885,     0.335974246, 
	    0.368491828,     0.404103488,     0.443095446,     0.485779136, 
	    0.532493234,     0.583605707,     0.639516354,     0.700658977, 
	    0.767504096,     0.840561509,     0.920383334,      1.00756681, 
	     1.10275745,      1.20665205,      1.32000256,      1.44361889, 
	     1.57837236,      1.72519934,      1.88510525,      2.05916667, 
	     2.24853587,      2.45444322,      2.67820048,      2.92120337, 
	     3.18493319,      3.47095847,      3.78093672,      4.11661386, 
	     4.47982359,      4.87248468,      5.29659796,      5.75424194, 
	     6.24756479,      6.77877522,      7.35013199,       7.9639287, 
	     8.62247753,       9.3280859,      10.0830355,      10.8895483, 
	     11.7497616,      12.6656818,      13.6391459,       14.671772, 
	     15.7649031,      16.9195518,      18.1363354,      19.4154034, 
	     20.7563629,      22.1581993,      23.6191978,      25.1368542, 
	     26.7077999,      28.3277054,      29.9912109,      31.6918564, 
	     33.4220123,      35.1728401,      36.9342651,      38.6949577, 
	     40.4423866,      42.1628456,      43.8415756,      45.4628983, 
	     47.0104179,      48.4672623,      49.8164101,      51.0410385, 
	     52.1249695,      53.0531502,      53.8122063,      54.3910103, 
	     54.7813034,      54.9783134,      54.9813919,      54.7945633, 
	     54.4270477,      53.8936539,      53.2150497,        52.41782, 
	     51.5343208,      50.6022682,      49.6640358,      48.7656441, 
	       47.95541,      47.2823181,      46.7940636,      46.5348778, 
	     46.5432014,       46.849247,      47.4726562,       48.420311, 
	     49.6845016,      51.2416191,      53.0514679,      55.0574455, 
	      57.187603,       59.356739,      61.4695396,       63.424675, 
	     65.1198349,      66.4573746,      67.3504486,       67.729126, 
	     67.5461884,      66.7820969,      65.4486618,      63.5909958, 
	     61.2873077,      58.6463623,      55.8024063,      52.9077072, 
	     50.1229668,      47.6062126,      45.5008659,       43.923996, 
	     42.9557343,      42.6309929,      42.9343758,      43.7990532, 
	     45.1100044,       46.711586,      48.4189796,      50.0326614, 
	     51.3545609,      52.2045021,      52.4352913,      51.9450874, 
	     50.6857491,      48.6664925,      45.9524689,      42.6585007, 
	     38.9386177,      34.9724808,      30.9500198,      27.0557327, 
	     23.4541149,      20.2775459,      17.6176262,      15.5207052, 
	     13.9877377,      12.9783249,       12.418232,       12.209404, 
	     12.2412558,       12.401947,      12.5884552,      12.7145233, 
	     12.7158709,      12.5524969,      12.2082186,      11.6879463, 
	     11.0133543,      10.2176933,      9.34045506,      8.42250156, 
	     7.50204515,      6.61176348,      5.77706861,      5.01546812, 
	      4.3368206,      3.74425268,      3.23547411,      2.80425739, 
	      2.4418931,      2.13847399,      1.88393486,      1.66880786, 
	     1.48471224,      1.32460427,      1.18284678,      1.05514359, 
	    0.938394129,     0.830503702,     0.730181158,     0.636742115, 
	    0.549929857,      0.46975863,     0.396382481,     0.329989672, 
	    0.270721495,     0.218615383,     0.173570439,     0.135333925, 
	    0.103505999,    0.0775592849,    0.0568695553,    0.0407526605, 
	   0.0285033733,    0.0194318648,     0.012894419,   0.00831620116, 
	  0.00520504266,   0.00315654068,   0.00185168535,   0.00104891404, 
	 0.000572720193,  0.000300850486,  0.000151740125,  7.33306879e-05, 
	 3.38809587e-05,  1.49318385e-05,  6.26197652e-06, 
/* Z=44 -> Ru */
	              0,    0.0323932245,    0.0355855227,    0.0390909836, 
	   0.0429401062,    0.0471663214,    0.0518062785,    0.0569001362, 
	   0.0624919087,    0.0686298013,    0.0753666312,    0.0827602297, 
	   0.0908739269,    0.0997770429,     0.109545447,     0.120262168, 
	    0.132018015,     0.144912317,     0.159053653,     0.174560741, 
	    0.191563249,     0.210202873,     0.230634317,     0.253026426, 
	    0.277563483,     0.304446429,     0.333894402,     0.366146207, 
	    0.401461869,     0.440124601,     0.482442468,     0.528750539, 
	    0.579412937,     0.634825289,     0.695416749,     0.761653125, 
	    0.834039032,     0.913120985,      0.99949044,      1.09378672, 
	      1.1967001,      1.30897558,      1.43141603,      1.56488526, 
	     1.71031272,      1.86869526,      2.04110241,       2.2286787, 
	     2.43264675,       2.6543119,      2.89506245,      3.15637445, 
	     3.43981218,      3.74702954,      4.07977104,      4.43987083, 
	     4.82925081,      5.24991894,      5.70396376,      6.19354868, 
	      6.7209034,      7.28831339,      7.89810658,       8.5526371, 
	     9.25426483,      10.0053349,      10.8081474,      11.6649246, 
	     12.5777798,      13.5486717,      14.5793581,      15.6713419, 
	     16.8258171,      18.0436058,      19.3250828,      20.6701088, 
	     22.0779476,      23.5471821,      25.0756416,      26.6602993, 
	     28.2971935,      29.9813519,      31.7067051,      33.4660263, 
	     35.2508698,      37.0515404,      38.8570747,      40.6552391, 
	     42.4326057,      44.1745872,       45.865612,       47.489254, 
	     49.0284882,      50.4659729,      51.7843895,      52.9668617, 
	     53.9974213,      54.8615341,      55.5466843,       56.042984, 
	     56.3438034,      56.4464149,      56.3526039,      56.0692215, 
	     55.6086502,       54.989151,      54.2350273,      53.3765793, 
	     52.4498367,      51.4959564,      50.5603752,      49.6915817, 
	     48.9395752,      48.3540497,      47.9822388,      47.8666077, 
	      48.042408,      48.5351982,      49.3585167,      50.5117912, 
	     51.9787407,      53.7263222,      55.7044868,      57.8467941, 
	     60.0720901,      62.2871857,      64.3906174,      66.2773819, 
	     67.8444901,      68.9970627,      69.6546707,      69.7575607, 
	     69.2721786,      68.1956482,      66.5586624,         64.4263, 
	     61.8965797,      59.0963554,      56.1747322,       53.294075, 
	     50.6191521,      48.3050385,      46.4847374,      45.2575722, 
	     44.6793861,      44.7557335,      45.4388351,      46.6289978, 
	     48.1805725,      49.9122963,      51.6211929,      53.0989227, 
	     54.1490631,       54.603775,      54.3381157,       53.280735, 
	     51.4198074,      48.8036995,      45.5363235,      41.7676888, 
	     37.6805649,      33.4746666,      29.3497791,      25.4895458, 
	     22.0473042,      19.1353474,      16.8184223,      15.1119699, 
	     13.9849873,      13.3670206,      13.1582546,      13.2414417, 
	     13.4942656,      13.8007746,      14.0607624,       14.196331, 
	     14.1553068,      13.9115753,      13.4628096,      12.8262644, 
	     12.0334358,      11.1244125,      10.1425362,      9.12991714, 
	     8.12407207,       7.1557765,      6.24806356,      5.41619921, 
	     4.66836977,      4.00683594,      3.42929935,      2.93028378, 
	     2.50239205,      2.13734078,      1.82674682,      1.56266344, 
	     1.33789825,      1.14616024,     0.982081831,      0.84116596, 
	    0.719690084,     0.614597917,     0.523392439,     0.444040537, 
	    0.374891102,     0.314605147,     0.262096673,     0.216480896, 
	    0.177028775,     0.143126637,     0.114241242,    0.0898911953, 
	   0.0696248636,    0.0530057512,    0.0396049395,    0.0289999172, 
	    0.020778507,     0.014546046,    0.0099339541,   0.00660796277, 
	  0.00427452894,   0.00268455781,   0.00163413538,  0.000962449063, 
	 0.000547465635,  0.000300200365,  0.000158377166,  8.02262657e-05, 
	 3.89367233e-05,  1.80657298e-05,  7.99455302e-06, 
/* Z=45 -> Rh */
	              0,    0.0355234817,    0.0390138514,    0.0428455621, 
	   0.0470517427,    0.0516686998,    0.0567362309,    0.0622979105, 
	   0.0684014931,    0.0750992596,    0.0824484825,    0.0905118212, 
	    0.099357903,       0.1090618,     0.119705647,     0.131379306, 
	    0.144181013,      0.15821816,     0.173608154,     0.190479234, 
	    0.208971485,     0.229237854,     0.251445323,     0.275775969, 
	    0.302428544,     0.331619561,     0.363585114,      0.39858222, 
	    0.436890781,     0.478815347,     0.524687111,     0.574866056, 
	    0.629743218,     0.689742982,     0.755325854,     0.826990843, 
	    0.905278563,     0.990774035,      1.08410978,      1.18596911, 
	     1.29708934,      1.41826534,      1.55035317,      1.69427335, 
	     1.85101438,      2.02163672,      2.20727634,      2.40914726, 
	     2.62854552,      2.86685205,      3.12553549,      3.40615296, 
	     3.71035337,      4.03987741,      4.39655733,      4.78231573, 
	     5.19916439,      5.64919949,      6.13459539,      6.65759897, 
	     7.22051859,      7.82571268,      8.47557259,      9.17250824, 
	     9.91892052,      10.7171803,      11.5695963,      12.4783792, 
	     13.4456043,      14.4731617,      15.5627136,      16.7156296, 
	     17.9329281,        19.21521,      20.5625858,      21.9745922, 
	     23.4501228,      24.9873257,      26.5835266,      28.2351456, 
	     29.9375992,      31.6852303,      33.4712257,      35.2875633, 
	      37.124958,      38.9728317,      40.8193207,      42.6512947, 
	     44.4544182,      46.2132683,      47.9114723,       49.531929, 
	     51.0570564,      52.4691238,      53.7506447,      54.8848228, 
	     55.8560524,      56.6505051,      57.2567177,      57.6662521, 
	     57.8743248,      57.8804741,      57.6891365,      57.3101883, 
	     56.7593422,       56.058445,      55.2355042,      54.3245621, 
	     53.3652382,      52.4019928,      51.4830551,      50.6590271, 
	     49.9811325,      49.4991913,      49.2593384,       49.301548, 
	     49.6571045,       50.346138,      51.3753166,      52.7359657, 
	     54.4026566,      56.3325653,      58.4656258,      60.7257233, 
	     63.0229263,      65.2568283,      67.3209381,      69.1080017, 
	      70.516037,      71.4548187,      71.8523636,      71.6610336, 
	     70.8627548,      69.4727783,      67.5416489,      65.1547852, 
	     62.4295692,       59.509716,      56.5570602,      53.7411613, 
	     51.2272453,      49.1633835,      47.6679535,        46.81847, 
	     46.6430016,      47.1151543,      48.1534424,      49.6254349, 
	     51.3565941,      53.1433487,      54.7692604,      56.0229568, 
	     56.7161293,      56.6998558,       55.877655,      54.2139587, 
	     51.7371597,      48.5369568,      44.7562256,      40.5783119, 
	     36.2109337,      31.8682976,      27.7532063,      24.0407715, 
	     20.8653164,      18.3116341,      16.4112759,      15.1440163, 
	     14.4441147,      14.2103968,      14.3188629,      14.6362629, 
	     15.0330954,      15.3946686,      15.6292496,       15.672761, 
	     15.4900188,      15.0728989,      14.4361467,      13.6117153, 
	     12.6425314,      11.5764751,       10.461175,      9.34000111, 
	     8.24939156,      7.21746683,      6.26376772,      5.39983559, 
	      4.6303606,      3.95461941,      3.36797071,      2.86324811, 
	     2.43193412,      2.06507683,      1.75393641,      1.49039626, 
	     1.26718044,      1.07793009,     0.917188108,     0.780329943, 
	    0.663472056,     0.563376248,       0.4773601,      0.40321669, 
	    0.339145005,     0.283686906,     0.235670149,     0.194155037, 
	     0.15838474,     0.127739415,     0.101695769,    0.0797933489, 
	   0.0616082773,    0.0467356332,    0.0347797982,    0.0253520161, 
	   0.0180734936,     0.012582045,   0.00854022522,    0.0056430758, 
	  0.00362407649,   0.00225840951,   0.00136332528,  0.000795851171, 
	 0.000448451086,   0.00024346383,  0.000127098465,  6.36716577e-05, 
	 3.05439098e-05,  1.39992799e-05,   6.1160772e-06, 
/* Z=46 -> Pd */
	              0,     0.038911555,    0.0427231491,    0.0469063073, 
	   0.0514969751,    0.0565345287,    0.0620621033,    0.0681269318, 
	   0.0747807473,     0.082080178,    0.0900872126,    0.0988696814, 
	    0.108501822,     0.119064808,     0.130647451,     0.143346846, 
	     0.15726909,     0.172530189,     0.189256817,     0.207587376, 
	    0.227672905,     0.249678284,     0.273783416,     0.300184488, 
	    0.329095334,     0.360749006,     0.395399392,     0.433322787, 
	    0.474819928,     0.520217836,     0.569871962,     0.624168396, 
	    0.683526397,     0.748400629,     0.819284081,     0.896710873, 
	    0.981258988,      1.07355344,       1.1742698,      1.28413701, 
	     1.40394127,      1.53452945,      1.67681229,      1.83176887, 
	     2.00044942,      2.18397951,       2.3835628,       2.6004858, 
	     2.83611894,      3.09192109,      3.36944151,      3.67032075, 
	     3.99629307,      4.34918642,      4.73092079,      5.14350843, 
	     5.58904886,      6.06972504,      6.58779573,      7.14558887, 
	     7.74548674,      8.38991547,      9.08132553,      9.82217312, 
	     10.6148911,      11.4618683,      12.3654079,      13.3276939, 
	     14.3507481,      15.4363766,      16.5861225,      17.8012009, 
	     19.0824261,      20.4301586,      21.8442097,      23.3237686, 
	     24.8673153,      26.4725342,      28.1362228,      29.8542004, 
	     31.6212368,       33.430954,      35.2757683,       37.146843, 
	     39.0340195,      40.9258385,      42.8095284,      44.6710625, 
	     46.4952316,      48.2657928,      49.9656448,      51.5770569, 
	     53.0819893,      54.4624367,      55.7008667,      56.7807045, 
	     57.6868935,      58.4064827,      58.9292755,      59.2484779, 
	     59.3613815,      59.2699814,      58.9815674,       58.509201, 
	     57.8720703,      57.0956726,      56.2117577,      55.2580414, 
	     54.2775993,      53.3179512,      52.4297867,      51.6653786, 
	     51.0766296,      50.7128716,      50.6184196,      50.8300133, 
	     51.3742142,      52.2649689,      53.5014267,      55.0662498, 
	      56.924511,      59.0234451,      61.2931175,      63.6481781, 
	     65.9907074,      68.2142258,      70.2086487,      71.8662186, 
	     73.0879364,      73.7902908,      73.9117661,      73.4187546, 
	      72.310173,      70.6204758,      68.4203949,      65.8152695, 
	     62.9405098,      59.9544144,      57.0284271,       54.335392, 
	     52.0365868,      50.2684631,      49.1302948,      48.6739235, 
	     48.8967361,      49.7388687,       51.085228,      52.7725334, 
	     54.6010284,      56.3500557,      57.7961731,      58.7321396, 
	     58.9850159,      58.4315453,      57.0092735,      54.7223129, 
	     51.6410904,      47.8961372,      43.6665154,      39.1640244, 
	     34.6147499,      30.2397461,      26.2367573,      22.7646999, 
	     19.9323502,      17.7922916,      16.3404198,      15.5209093, 
	      15.235734,      15.3575077,      15.7439747,      16.2524395, 
	     16.7525692,      17.1363297,      17.3243046,      17.2682152, 
	     16.9499168,      16.3776016,      15.5800924,      14.6002502, 
	     13.4883652,      12.2962904,      11.0727329,      9.85997772, 
	     8.69199753,      7.59381628,      6.58184862,      5.66491699, 
	     4.84564734,      4.12199545,      3.48870969,      2.93859935, 
	     2.46355581,      2.05529737,      1.70586574,      1.40791595, 
	     1.15484226,     0.940799356,     0.760653436,     0.609903514, 
	    0.484592825,      0.38122642,     0.296702176,     0.228257239, 
	    0.173428267,     0.130023301,    0.0961013734,    0.0699568987, 
	   0.0501061864,    0.0352743752,    0.0243814848,    0.0165270865, 
	   0.0109735159,    0.0071278899,    0.0045234249,   0.00280067488, 
	  0.00168933335,  0.000991205452,  0.000564825314,  0.000312061398, 
	  0.00016687077,  8.62060115e-05,  4.29413885e-05,  2.05835622e-05, 
	  9.4744637e-06,  4.17847286e-06,  1.76157675e-06,  7.08198684e-07, 
	 2.70814212e-07,  9.82423387e-08,  3.37153949e-08, 
/* Z=47 -> Ag */
	              0,    0.0425906554,    0.0467495658,    0.0513125695, 
	   0.0563186295,    0.0618104078,    0.0678346083,    0.0744423792, 
	   0.0816897079,    0.0896378532,    0.0983538553,      0.10791105, 
	    0.118389659,     0.129877388,     0.142470136,     0.156272694, 
	    0.171399593,     0.187975898,     0.206138223,     0.226035669, 
	    0.247831002,     0.271701723,      0.29784146,     0.326461256, 
	    0.357791096,      0.39208141,     0.429604888,     0.470658273, 
	    0.515564203,     0.564673424,     0.618367016,     0.677058578, 
	     0.74119699,     0.811268806,      0.88780117,     0.971364856, 
	     1.06257701,      1.16210485,      1.27066875,      1.38904524, 
	     1.51807153,      1.65864837,      1.81174409,      1.97839797, 
	     2.15972447,      2.35691619,      2.57124758,      2.80407858, 
	     3.05685687,      3.33112121,      3.62850285,      3.95072794, 
	     4.29961729,      4.67708635,        5.085145,      5.52589369, 
	     6.00152016,      6.51429272,      7.06655359,      7.66070652, 
	     8.29920673,      8.98454285,      9.71921635,       10.505723, 
	     11.3465223,      12.2440033,      13.2004547,      14.2180223, 
	     15.2986517,      16.4440517,      17.6556206,      18.9343967, 
	     20.2809696,      21.6954269,      23.1772499,      24.7252464, 
	     26.3374577,      28.0110607,      29.7422905,      31.5263386, 
	     33.3572807,      35.2279968,      37.1300964,       39.053894, 
	     40.9883575,      42.9211197,      44.8384972,      46.7255783, 
	     48.5662994,      50.3436356,      52.0397987,      53.6365204, 
	     55.1153755,      56.4581985,      57.6475487,      58.6672325, 
	      59.502903,      60.1426888,      60.5778542,      60.8034821, 
	     60.8191376,      60.6294975,      60.2448921,      59.6817474, 
	     58.9628525,      58.1174355,       57.181015,      56.1949196, 
	     55.2055397,      54.2632103,      53.4207573,      52.7316628, 
	     52.2479858,      52.0179253,      52.0832901,      52.4768181, 
	     53.2195969,      54.3186455,      55.7649231,      57.5318642, 
	     59.5746574,      61.8304253,      64.2194519,      66.6475296, 
	     69.0094681,      71.1937103,      73.0879593,      74.5855484, 
	     75.5922546,      76.0332031,      75.8592834,      75.0526428, 
	     73.6306915,      71.6480408,      69.1960907,      66.3997192, 
	      63.411274,      60.4016762,      57.5492287,      55.0266685, 
	     52.9874496,      51.5523529,      50.7976418,      50.7460785, 
	     51.3618393,      52.5502548,      54.1627121,      56.0067329, 
	     57.8605461,       59.491024,      60.6734161,      61.2110252, 
	     60.9528847,      59.8076477,      57.7522736,      54.8345451, 
	     51.1691513,      46.9276505,      42.3233147,      37.5923004, 
	     32.9730949,      28.6861286,      24.9156914,      21.7958355, 
	     19.4015789,      17.7462082,      16.7846375,      16.4222603, 
	     16.5279579,      16.9496346,       17.530302,      18.1229286, 
	     18.6024647,      18.8740768,      18.8771267,      18.5850906, 
	     18.0020733,      17.1568851,      16.0958118,      14.8751221, 
	     13.5541773,      12.1897821,      10.8320684,      9.52198696, 
	     8.29025459,      7.15750694,      6.13529062,      5.22757673, 
	     4.43248463,      3.74398804,      3.15343571,      2.65080309, 
	     2.22564268,       1.8677454,      1.56755853,      1.31640971, 
	     1.10659516,     0.931379676,     0.784947336,     0.662328541, 
	    0.559320867,     0.472410411,     0.398697615,     0.335825741, 
	    0.281911105,     0.235474363,     0.195371807,     0.160728425, 
	    0.130874246,     0.105286352,    0.0835390687,    0.0652638376, 
	    0.050119888,    0.0377756916,    0.0279001966,    0.0201620497, 
	   0.0142344711,   0.00980335381,   0.00657624472,   0.00429026829, 
	  0.00271776784,   0.00166902004,  0.000991991372,  0.000569640368, 
	 0.000315470272,   0.00016817663,  8.61329681e-05,  4.22938392e-05, 
	 1.98679409e-05,  8.90871434e-06,  3.80393953e-06, 
/* Z=48 -> Cd */
	              0,    0.0465783365,    0.0511120185,    0.0560847335, 
	    0.061538659,    0.0675199553,    0.0740791708,    0.0812716037, 
	    0.089157775,    0.0978038907,     0.107282378,     0.117672421, 
	    0.129060626,       0.1415416,      0.15521881,     0.170205235, 
	    0.186624303,     0.204610735,      0.22431168,     0.245887607, 
	    0.269513607,     0.295380652,     0.323696792,     0.354688823, 
	    0.388603687,     0.425710291,     0.466301113,     0.510694325, 
	    0.559235752,     0.612301111,     0.670298338,     0.733669937, 
	    0.802895844,     0.878495991,     0.961033404,      1.05111694, 
	     1.14940524,      1.25660884,      1.37349498,      1.50089025, 
	      1.6396842,      1.79083395,      1.95536721,      2.13438582, 
	     2.32907009,      2.54068208,      2.77056885,      3.02016616, 
	     3.29100013,      3.58469129,      3.90295482,      4.24760246, 
	       4.620543,       5.0237813,      5.45941544,      5.92963409, 
	     6.43671179,      6.98300171,      7.57092476,      8.20295906, 
	     8.88162708,      9.60947227,      10.3890448,       11.222868, 
	     12.1134138,      13.0630693,        14.07409,      15.1485624, 
	     16.2883472,      17.4950275,      18.7698441,      20.1136246, 
	     21.5267143,      23.0088921,      24.5592861,      26.1762886, 
	     27.8574657,      29.5994568,      31.3978977,      33.2473145, 
	     35.1410675,      37.0712471,      39.0286484,      41.0027046, 
	     42.9814987,      44.9517441,      46.8988647,      48.8070641, 
	     50.6594658,      52.4383011,      54.1251678,      55.7013245, 
	     57.1480827,      58.4472351,      59.5815697,      60.5354424, 
	     61.2954063,      61.8508492,      62.1947174,      62.3241692, 
	     62.2412643,      61.9535484,       61.474575,      60.8242683, 
	     60.0291214,       59.122139,      58.1425629,      57.1352272, 
	     56.1496353,      55.2386398,      54.4567757,      53.8582573, 
	     53.4946823,      53.4124603,      53.6501541,      54.2357597, 
	     55.1841278,      56.4946785,      58.1495781,      60.1126137, 
	     62.3288422,       64.725296,      67.2127991,      69.6889114, 
	     72.0421219,      74.1571121,      75.9209137,      77.2297668, 
	     77.9961853,      78.1559067,      77.6740799,      76.5502396, 
	      74.821434,      72.5630646,      69.8869858,      66.9367142, 
	     63.8796692,      60.8967094,      58.1695442,      55.8667984, 
	     54.1298103,      53.0594292,      52.7051048,      53.0575256, 
	      54.045887,      55.5404129,       57.360363,      59.2871704, 
	     61.0817146,      62.5043716,      63.3359375,      63.3973885, 
	     62.5665627,      60.7899361,      58.0882835,      54.5554924, 
	     50.3506546,      45.6840744,      40.7985954,      35.9480743, 
	     31.3751621,      27.2906094,      23.8561268,      21.1725845, 
	     19.2745514,      18.1315994,      17.6559601,      17.7153683, 
	     18.1493855,      18.7871113,      19.4641724,      20.0371361, 
	     20.3939571,      20.4597931,      20.1981525,      19.6079273, 
	     18.7173901,      17.5763092,      16.2474613,      14.7985926, 
	     13.2955761,      11.7973061,      10.3523884,      8.99758053, 
	     7.75768948,      6.64654922,      5.66870117,      4.82140303, 
	     4.09668541,      3.48323488,      2.96799064,      2.53738952, 
	       2.178267,      1.87845182,       1.6271069,      1.41488171, 
	     1.23393011,      1.07784128,     0.941518128,     0.821025491, 
	    0.713426292,     0.616612792,     0.529140592,      0.45007056, 
	     0.37882176,     0.315041035,     0.258491933,     0.208966076, 
	    0.166218713,      0.12992847,    0.0996793956,    0.0749626681, 
	   0.0551930778,    0.0397356115,    0.0279367846,    0.0191560257, 
	   0.0127932243,   0.00830983091,   0.00524215167,   0.00320680556, 
	  0.00189929362,   0.00108729268,  0.000600591826,  0.000319519721, 
	 0.000163406017,  8.01707938e-05,  3.76550342e-05,  1.68936513e-05, 
	 7.22276945e-06,  2.93561288e-06,  1.13133365e-06, 
/* Z=49 -> In */
	              0,    0.0508931205,    0.0558304451,    0.0612442531, 
	   0.0671801493,    0.0736880228,    0.0808224976,     0.088643305, 
	    0.097215794,     0.106611431,     0.116908342,     0.128191948, 
	    0.140555605,     0.154101297,      0.16894047,     0.185194775, 
	    0.202997074,     0.222492352,     0.243838802,     0.267208964, 
	    0.292791009,     0.320789963,     0.351429254,     0.384952217, 
	    0.421623677,     0.461731881,       0.5055902,     0.553539455, 
	    0.605949581,     0.663222611,     0.725794554,     0.794138253, 
	    0.868766069,     0.950232923,      1.03913915,      1.13613355, 
	     1.24191737,      1.35724711,      1.48293805,      1.61986864, 
	     1.76898324,      1.93129694,      2.10789847,      2.29995441, 
	     2.50871301,      2.73550749,      2.98175955,      3.24898171, 
	     3.53878069,      3.85285974,      4.19301939,       4.5611577, 
	     4.95927238,       5.3894558,      5.85389566,      6.35486746, 
	     6.89473057,      7.47591972,      8.10093212,      8.77231598, 
	     9.49265385,      10.2645407,      11.0905619,      11.9732628, 
	     12.9151192,      13.9184952,      14.9856033,      16.1184559, 
	     17.3188038,      18.5880909,      19.9273701,      21.3372402, 
	     22.8177681,      24.3684006,      25.9878788,      27.6741505, 
	     29.4242687,      31.2343102,      33.0992699,      35.0129929, 
	     36.9680672,      38.9557991,      40.9661255,      42.9876022, 
	     45.0074005,      47.0113335,      48.9839325,      50.9085503, 
	     52.7675323,      54.5424347,      56.2143097,      57.7640457, 
	     59.1727982,      60.4224548,      61.4961929,       62.379097, 
	     63.0587921,      63.5261497,      63.7759895,      63.8077621, 
	     63.6262093,      63.2419395,      62.6718788,      61.9395599, 
	     61.0752335,      60.1156845,      59.1037865,      58.0877113, 
	     57.1198044,      56.2550583,       55.549263,      55.0567856, 
	     54.8281097,      54.9071426,      55.3284416,      56.1144981, 
	     57.2732048,      58.7957306,      60.6549568,      62.8046722, 
	     65.1797333,      67.6972809,      70.2591553,      72.7555542, 
	     75.0698242,      77.0843735,      78.6873398,      79.7797928, 
	     80.2830429,      80.1454086,      79.3480835,      77.9094009, 
	     75.8869247,      73.3770218,      70.5114136,      67.4507065, 
	     64.3749771,      61.4717865,      58.9223328,      56.8867722, 
	     55.4898834,      54.8084183,      54.8615227,      55.6054535, 
	     56.9334717,      58.6814461,      60.6390419,      62.5658417, 
	     64.2111282,      65.3356171,      65.7329712,      65.2491226, 
	     63.7972145,      61.3666458,      58.0251579,      53.9135818, 
	     49.2337036,      44.2303505,      39.1694832,      34.3144646, 
	     29.9028778,      26.1263428,      23.1152878,      20.9302788, 
	     19.5606384,      18.9302464,      18.9096394,      19.3326817, 
	     20.0156631,      20.7764301,      21.4513741,      21.9084892, 
	     22.0554848,      21.8426323,      21.2607861,      20.3355656, 
	     19.1190128,      17.6801434,      16.0956116,      14.4415693, 
	     12.7872448,      11.1906433,      9.69622231,      8.33437347, 
	     7.12227869,      6.06569481,      5.16125917,      4.39894533, 
	     3.76439285,       3.2409389,      2.81125712,      2.45858097, 
	     2.16753888,      1.92465425,      1.71857929,      1.54012728, 
	     1.38216376,      1.23940492,      1.10816097,     0.986054242, 
	    0.871733069,      0.76459831,     0.664555252,     0.571799517, 
	    0.486642241,     0.409377515,     0.340191364,     0.279108435, 
	    0.225971714,     0.180447757,      0.14205049,      0.11017625, 
	   0.0841434449,    0.0632320791,     0.046719145,    0.0339078307, 
	   0.0241495185,    0.0168586355,    0.0115209911,   0.00769672729, 
	  0.00501899654,   0.00318948762,   0.00197184901,   0.00118383672, 
	 0.000688897911,  0.000387801381,  0.000210750528,   0.00011033527, 
	 5.55258812e-05,   2.6799542e-05,  1.23762711e-05, 
/* Z=50 -> Sn */
	              0,    0.0555638708,    0.0609360971,    0.0668249503, 
	   0.0732796937,    0.0803542063,    0.0881074294,    0.0966038257, 
	     0.10591387,     0.116114631,     0.127290294,     0.139532909, 
	       0.152943,      0.16763036,     0.183714896,     0.201327473, 
	    0.220610932,     0.241721079,     0.264827877,     0.290116578, 
	    0.317789108,     0.348065466,     0.381185174,     0.417409062, 
	    0.457020789,     0.500328958,     0.547668934,      0.59940511, 
	    0.655933142,     0.717682302,     0.785118222,     0.858745396, 
	    0.939110279,      1.02680409,      1.12246621,      1.22678733, 
	     1.34051275,      1.46444666,      1.59945464,      1.74646866, 
	     1.90648985,      2.08059311,       2.2699306,      2.47573543, 
	      2.6993258,      2.94210744,      3.20557833,      3.49133039, 
	     3.80105209,       4.1365304,      4.49965191,      4.89240265, 
	     5.31686783,      5.77522945,      6.26976204,      6.80282784, 
	     7.37686872,      7.99439764,      8.65798473,      9.37024403, 
	      10.133812,      10.9513273,      11.8254013,      12.7585917, 
	     13.7533646,      14.8120508,       15.936801,      17.1295357, 
	     18.3918839,      19.7251148,      21.1300793,      22.6071186, 
	     24.1559925,      25.7757893,      27.4648285,      29.2205811, 
	     31.0395603,      32.9172363,      34.8479462,      36.8247986, 
	     38.8396187,      40.8828697,      42.9436264,      45.0095482, 
	     47.0669022,      49.1006012,      51.0943031,      53.0305405, 
	     54.8909264,      56.6564102,      58.3075752,      59.8250504, 
	      61.189949,      62.3844032,      63.3921471,      64.1991501, 
	     64.7943192,      65.1702042,      65.3237076,      65.2567673, 
	     64.9770126,      64.4982681,      63.8409729,      63.0323677, 
	     62.1064796,      61.1038399,      60.0708427,      59.0587883, 
	     58.1225586,      57.3188858,      56.7043037,      56.3327446, 
	     56.2529221,      56.5055008,      57.1202965,      58.1135597, 
	     59.4855423,      61.2185783,      63.2758255,      65.6008759, 
	     68.1183853,      70.7359238,      73.3470001,      75.8353958, 
	     78.0806046,      79.9643402,      81.3776398,      82.2284088, 
	     82.4486771,      82.0012283,       80.884903,      79.1379852, 
	     76.8391495,      74.1055298,      71.0877228,      67.9615631, 
	      64.917038,      62.1449013,      59.8217888,      58.0950661, 
	     57.0686722,      56.7914658,      57.2493629,      58.3624992, 
	     59.9880753,      61.9291916,      63.9492683,       65.791008, 
	     67.1983795,      67.9394913,      67.8282547,      66.7424698, 
	     64.6364212,      61.5466003,      57.5896873,      52.9529381, 
	     47.8777237,      42.6378326,      37.5147057,      32.7720947, 
	     28.6327763,      25.2597828,      22.7440052,      21.0995693, 
	     20.2671261,      20.1245461,      20.5034351,      21.2092667, 
	     22.0425739,      22.8186169,      23.3833961,      23.6244888, 
	     23.4760876,      22.9184437,      21.9725723,      20.6916485, 
	     19.1505985,      17.4354305,      15.6334934,      13.8255987, 
	     12.0803881,       10.451067,      8.97429276,      7.67081165, 
	     6.54735184,      5.59928036,      4.81355906,      4.17166328, 
	     3.65218902,      3.23301911,      2.89297438,      2.61295843, 
	      2.3766489,       2.1708045,      1.98527586,      1.81279635, 
	     1.64862955,      1.49013412,      1.33629715,      1.18727505, 
	     1.04397225,     0.907678723,     0.779776275,     0.661519349, 
	    0.553888142,     0.457506418,     0.372614473,     0.299084038, 
	    0.236463308,     0.184039906,     0.140911639,     0.106057331, 
	   0.0784016997,     0.056870956,    0.0404372253,    0.0281513613, 
	   0.0191647373,    0.0127411308,   0.00826019887,   0.00521420268, 
	   0.0031996374,    0.0019054357,   0.00109925715,  0.000613214972, 
	 0.000330139301,  0.000171190884,  8.53207748e-05,    4.078247e-05, 
	 1.86530888e-05,  8.14437863e-06,   3.3862716e-06, 
/* Z=51 -> Sb */
	              0,    0.0606186949,    0.0664593354,    0.0728596225, 
	   0.0798727125,    0.0875567272,    0.0959752202,     0.105197683, 
	    0.115300037,     0.126365304,     0.138484135,     0.151755556, 
	     0.16628772,     0.182198718,     0.199617416,     0.218684435, 
	    0.239553168,     0.262390912,      0.28737995,     0.314719081, 
	    0.344624728,     0.377332628,     0.413099349,     0.452204049, 
	    0.494950265,     0.541667998,     0.592715681,     0.648482621, 
	    0.709391177,      0.77589941,     0.848503828,     0.927742124, 
	     1.01419616,      1.10849535,      1.21131957,      1.32340312, 
	     1.44553745,      1.57857573,      1.72343624,      1.88110602, 
	     2.05264497,       2.2391901,      2.44195867,      2.66225195, 
	     2.90146017,      3.16106391,       3.4426384,      3.74785662, 
	     4.07848978,      4.43641043,      4.82359219,      5.24210835, 
	      5.6941328,      6.18193579,      6.70787811,      7.27440643, 
	     7.88404322,      8.53937721,      9.24304581,      9.99772263, 
	     10.8060913,      11.6708269,      12.5945597,      13.5798473, 
	     14.6291304,      15.7446928,      16.9286098,      18.1826897, 
	     19.5084133,      20.9068642,      22.3786564,      23.9238529, 
	     25.5418758,      27.2314186,      28.9903603,      30.8156528, 
	     32.7032394,       34.647953,      36.6434364,       38.682045, 
	     40.7547989,      42.8513107,      44.9597702,      47.0669403, 
	     49.1581726,      51.2174911,      53.2277107,      55.1705894, 
	     57.0270576,      58.7775192,       60.402195,      61.8815536, 
	     63.1968079,      64.3304749,      65.2670135,      65.9934692, 
	     66.5002213,      66.7817001,      66.8370895,      66.6710205, 
	     66.2941742,      65.7237625,      64.9838333,      64.1054001, 
	     63.1262932,      62.0906982,      61.0483818,      60.0535202, 
	     59.1631851,      58.4354134,      57.9269447,      57.6906548, 
	      57.772789,      58.2100601,      59.0268059,      60.2323112, 
	     61.8185654,      63.7585678,      66.0054092,       68.492363, 
	     71.1340485,      73.8289185,      76.4629288,      78.9145355, 
	     81.0608139,      82.7844391,      83.9812317,      84.5678406, 
	     84.4889603,      83.7235794,      82.2895584,      80.2460175, 
	     77.6929398,      74.7677155,      71.6383743,      68.4936676, 
	     65.5304184,      62.9388084,      60.8867264,      59.5044365, 
	     58.8710175,      59.0040588,      59.8539581,      61.3038254, 
	     63.1755562,      65.2420197,      67.2445984,      68.9148331, 
	     69.9981155,      70.2772903,      69.5936661,      67.8633118, 
	     65.0866013,      61.3500519,       56.819828,      51.7275124, 
	     46.3493729,      40.9811745,      35.9110985,      31.3935642, 
	     27.6267414,      24.7361794,      22.7662277,      21.6801605, 
	     21.3686733,      21.6655769,      22.3684654,      23.2617092, 
	     24.1388454,      24.8217983,      25.1749306,      25.1128292, 
	     24.6017227,      23.6551876,      22.3255901,      20.6928654, 
	     18.8524303,      16.9037437,      14.9405832,      13.0438557, 
	     11.2770615,      9.68434906,      8.29079914,      7.10435343, 
	     6.11884212,      5.31755209,      4.67685556,      4.16956282, 
	     3.76775956,      3.44501209,      3.17790914,       2.9469831, 
	     2.73709369,      2.53738117,      2.34090209,      2.14405632, 
	     1.94590414,      1.74745095,       1.5509665,      1.35937977, 
	      1.1757791,      1.00303054,      0.84351486,     0.698974133, 
	    0.570450664,     0.458301902,     0.362268656,     0.281578183, 
	    0.215065077,     0.161295876,     0.118686892,    0.0856079236, 
	     0.06046772,    0.0417790934,    0.0282038543,    0.0185789578, 
	   0.0119263874,   0.00744980341,   0.00452133035,   0.00266178953, 
	  0.00151750981,  0.000836311607,  0.000444707752,  0.000227722267, 
	 0.000112066482,  5.28883793e-05,  2.38831653e-05,  1.02957802e-05, 
	 4.22673111e-06,  1.64824439e-06,  6.08909602e-07, 
/* Z=52 -> Te */
	              0,    0.0660882443,    0.0724332556,    0.0793839917, 
	   0.0869977549,    0.0953371897,      0.10447076,     0.114473291, 
	    0.125426516,     0.137419745,     0.150550485,     0.164925218, 
	    0.180660173,     0.197882205,     0.216729671,     0.237353519, 
	    0.259918332,     0.284603447,     0.311604381,     0.341134012, 
	    0.373424113,      0.40872705,     0.447317302,      0.48949334, 
	    0.535579622,     0.585928619,     0.640922964,     0.700978041, 
	    0.766544342,      0.83810991,     0.916203678,      1.00139785, 
	     1.09431159,      1.19561374,       1.3060267,      1.42632973, 
	     1.55736244,      1.70002913,      1.85530221,      2.02422595, 
	     2.20792127,      2.40758824,      2.62451172,      2.86006355, 
	     3.11570692,      3.39299917,      3.69359517,      4.01924944, 
	     4.37181759,      4.75325775,      5.16563129,      5.61110067, 
	     6.09192657,      6.61046696,        7.169168,      7.77055788, 
	     8.41723824,        9.111866,      9.85714436,      10.6557999, 
	     11.5105572,      12.4241161,      13.3991175,      14.4381027, 
	     15.5434809,      16.7174702,      17.9620495,      19.2788982, 
	     20.6693287,      22.1342125,      23.6739044,      25.2881565, 
	     26.9760284,      28.7357922,      30.5648422,      32.4595909, 
	     34.4153786,      36.4263649,      38.4854698,      40.5842667, 
	     42.7129478,      44.8602638,      47.0135002,      49.1585121, 
	     51.2797585,      53.3603897,      55.3823967,      57.3268089, 
	     59.1739464,      60.9037514,        62.49617,      63.9316368, 
	      65.191597,      66.2591171,      67.1195374,      67.7611923, 
	     68.1761093,      68.3608017,      68.3169479,      68.0520401, 
	     67.5799942,      66.9215164,      66.1043777,      65.1634064, 
	     64.1401749,      63.0824394,      62.0431328,      61.0790176, 
	     60.2489738,      59.6118507,      59.2240219,      59.1366653, 
	     59.3928604,      60.0246315,      61.0501022,      62.4709511, 
	     64.2703171,      66.4114151,       68.837059,       71.470192, 
	     74.2157135,      76.9635239,      79.5929565,      81.9784164, 
	     83.9960709,      85.5313187,      86.4866714,      86.7894821, 
	      86.398941,       85.311821,      83.5661697,      81.2425308, 
	     78.4620667,      75.3814545,      72.1843262,      69.0696793, 
	     66.2377167,      63.8740616,      62.1335678,      61.1251678, 
	     60.8993149,      61.4394341,      62.6587372,      64.4032288, 
	     66.4611435,      68.5785217,      70.4797745,      71.8916321, 
	     72.5681152,      72.3141937,      71.0055466,      68.6022873, 
	     65.1548843,      60.8015747,      55.7570648,      50.2935982, 
	     44.7161331,      39.3341026,      34.4327087,      30.2468319, 
	     26.9403534,      24.5932121,      23.1975479,      22.6632214, 
	     22.8318367,      23.4972992,      24.4301453,      25.4025459, 
	     26.2108822,      26.6934204,      26.7414398,      26.3032303, 
	     25.3814259,      24.0248775,      22.3169022,      20.3618336, 
	     18.2716808,      16.1544132,      14.1047802,      12.1982059, 
	     10.4877462,      9.00378132,      7.75590181,      6.73630905, 
	     5.92407131,      5.28961992,      4.79901361,      4.41761923, 
	     4.11300516,      3.85697699,      3.62676859,      3.40549421, 
	     3.18199492,      2.95024514,      2.70847321,      2.45814705, 
	     2.20294118,       1.9477793,      1.69801366,      1.45877814, 
	     1.23452377,      1.02872968,     0.843768835,     0.680899918, 
	    0.540353894,     0.421483338,     0.322947204,     0.242907107, 
	     0.17921716,     0.129594266,    0.0917608291,    0.0635556728, 
	   0.0430129096,    0.0284102559,    0.0182906762,    0.0114623057, 
	  0.00698192045,    0.0041273823,    0.0023641251,   0.00130986481, 
	 0.000700757781,  0.000361313549,    0.0001791947,  8.53100501e-05, 
	 3.89024608e-05,  1.69542964e-05,  7.04509512e-06,  2.78438301e-06, 
	 1.04396725e-06,  3.70330866e-07,  1.23940467e-07, 
/* Z=53 -> I */
	              0,    0.0720055401,    0.0788934827,    0.0864365026, 
	   0.0946962982,      0.10374032,      0.11364226,     0.124482617, 
	    0.136349291,      0.14933826,     0.163554341,     0.179111883, 
	    0.196135685,     0.214761853,     0.235138848,     0.257428497, 
	    0.281807214,      0.30846712,     0.337617517,     0.369486332, 
	    0.404321462,     0.442392737,       0.4839935,     0.529442608, 
	    0.579086304,     0.633300841,     0.692494273,     0.757109165, 
	    0.827625394,     0.904562652,     0.988483548,      1.07999635, 
	     1.17975914,      1.28848195,      1.40693104,       1.5359329, 
	     1.67637682,      1.82922006,      1.99549115,      2.17629361, 
	     2.37281084,      2.58630872,      2.81814098,       3.0697515, 
	     3.34267879,      3.63855863,      3.95912623,       4.3062191, 
	     4.68177891,      5.08785057,      5.52658224,      6.00022507, 
	     6.51112747,      7.06173325,      7.65457201,      8.29225159, 
	     8.97744656,        9.712883,      10.5013208,      11.3455315, 
	     12.2482748,       13.212265,      14.2401419,      15.3344221, 
	     16.4974613,      17.7313995,      19.0381031,      20.4190998, 
	     21.8755093,      23.4079704,      25.0165482,      26.7006607, 
	     28.4589729,      30.2893047,      32.1885338,       34.152504, 
	     36.1759109,      38.2522278,      40.3736076,      42.5308304, 
	     44.7132301,      46.9086723,        49.10355,      51.2828064, 
	     53.4300041,      55.5274582,      57.5563736,      59.4971046, 
	     61.3294182,      63.0328903,      64.5873032,      65.9731827, 
	     67.1723557,      68.1685944,      68.9483109,      69.5012894, 
	     69.8214417,      69.9075317,      69.7639236,      69.4011688, 
	     68.8365479,      68.0943909,      67.2062531,      66.2107697, 
	     65.1532135,      64.0847397,      63.0612259,      62.1416931, 
	     61.3863907,      60.8544617,      60.6013222,      60.6757736, 
	     61.1170197,      61.9516716,      63.1909561,      64.8283005, 
	     66.8374939,      69.1716232,      71.7630768,      74.5246353, 
	     77.3518524,      80.1268616,      82.7234421,      85.0132828, 
	      86.873291,      88.1934357,      88.8848343,      88.8874359, 
	     88.1766739,      86.7684937,      84.7220535,      82.1395569, 
	     79.1628647,      75.9665604,      72.7476501,      69.7123489, 
	     67.0605011,      64.9690094,      63.5753975,      62.9632492, 
	     63.1509781,      64.0855179,      65.6419983,       67.630043, 
	     69.8067703,      71.8956757,      73.6100082,      74.6785507, 
	     74.8713989,      74.0229721,      72.0497894,      68.9609528, 
	     64.8598785,      59.9367828,      54.4524002,      48.7144165, 
	     43.0488358,      37.7693329,       33.147789,      29.3892918, 
	     26.6143475,      24.8504219,      24.0335674,      24.0198631, 
	     24.6049843,      25.5492649,      26.6049404,      27.5421143, 
	     28.1703815,      28.3538532,      28.0184765,      27.1516209, 
	     25.7950497,      24.0329895,      21.9774799,      19.7531414, 
	     17.4831219,      15.2776909,      13.2261333,       11.392231, 
	     9.81311989,      8.50094128,      7.44656944,      6.62461281, 
	      5.9989109,      5.52788639,      5.16925669,         4.88378, 
	     4.63789606,      4.40525198,      4.16722345,      3.91261888, 
	     3.63678336,      3.34033442,       3.0277307,      2.70585585, 
	     2.38273501,      2.06647158,      1.76444387,      1.48276627, 
	     1.22599244,     0.997026861,     0.797194123,     0.626423657, 
	    0.483500898,     0.366350144,     0.272316933,     0.198428646, 
	    0.141618535,    0.0989050344,    0.0675239488,    0.0450151786, 
	   0.0292688124,    0.0185371339,    0.0114203822,   0.00683426904, 
	  0.00396655919,   0.00222918531,   0.00121101737,  0.000634824915, 
	 0.000320509833,  0.000155547139,  7.24139027e-05,  3.22689266e-05, 
	 1.37330881e-05,   5.5685814e-06,   2.1460321e-06,  7.83998132e-07, 
	 2.70769249e-07,  8.81563054e-08,  2.69764371e-08, 
/* Z=54 -> Xe */
	              0,    0.0784064084,    0.0858786702,    0.0940588042, 
	     0.10301324,     0.112814508,     0.123541869,     0.135281831, 
	    0.148128852,     0.162186012,     0.177565813,     0.194390997, 
	    0.212795407,     0.232924968,     0.254938811,     0.279010266, 
	     0.30532819,     0.334098279,     0.365544409,     0.399910212, 
	    0.437460691,     0.478484035,     0.523293436,     0.572229028, 
	    0.625660241,     0.683987975,     0.747646928,     0.817108393, 
	    0.892882884,     0.975522995,      1.06562662,      1.16383994, 
	     1.27086079,      1.38744271,      1.51439774,      1.65260065, 
	     1.80299306,      1.96658635,      2.14446759,      2.33780098, 
	      2.5478344,      2.77590179,      3.02342701,      3.29192853, 
	     3.58302116,      3.89842033,       4.2399435,      4.60951328, 
	     5.00915575,      5.44100428,      5.90729475,      6.41036558, 
	     6.95265293,      7.53668404,      8.16506958,      8.84049511, 
	     9.56570435,      10.3434849,      11.1766491,      12.0680075, 
	     13.0203419,       14.036377,      15.1187344,      16.2698936, 
	     17.4921436,      18.7875252,      20.1577759,       21.604248, 
	     23.1278496,      24.7289543,      26.4073219,      28.1619968, 
	     29.9912224,      31.8923378,      33.8616753,      35.8944702, 
	     37.9847527,      40.1252632,      42.3073883,      44.5210686, 
	     46.7547684,      48.9954681,      51.2286415,      53.4383392, 
	     55.6072578,       57.716877,      59.7476807,      61.6793976, 
	     63.4913368,      65.1627884,      66.6735001,      68.0042038, 
	     69.1372757,      70.0573502,      70.7521057,      71.2129669, 
	     71.4359055,      71.4221649,      71.1789551,      70.7200394, 
	     70.0662155,      69.2455292,      68.2933502,      67.2520981, 
	     66.1706467,      65.1033783,      64.1088104,      63.2478714, 
	     62.5817108,        62.16922,      62.0642128,      62.3124428, 
	     62.9485207,      63.9929352,      65.4493561,      67.3023758, 
	     69.5159836,      72.0329208,      74.7751236,      77.6454239, 
	     80.5306549,      83.3060303,      85.8409653,      88.0059509, 
	     89.6803284,      90.7605133,      91.1681595,      90.8576355, 
	     89.8222122,      88.0981903,      85.7664566,      82.9508209, 
	     79.8129196,      76.5435028,      73.3503571,      70.4435043, 
	     68.0185165,      66.2393112,      65.2218857,      65.0206757, 
	     65.6191406,      66.9260406,      68.7782593,      70.9506989, 
	       73.17276,      75.1504898,      76.5923691,      77.2365646, 
	     76.8767319,      75.3837585,      72.7208405,      68.9500046, 
	      64.229126,      58.7991486,      52.9627686,      47.0563354, 
	     41.4180222,      36.3554802,      32.1166115,      28.8668022, 
	     26.6751442,      25.5113792,      25.2537003,      25.7063427, 
	     26.6244984,      27.7432671,      28.8067875,      29.5940056, 
	     29.9381275,      29.7379456,       28.960556,      27.6360989, 
	     25.8462563,      23.7087364,      21.3601532,      18.9395733, 
	     16.5744362,      14.3700781,      12.4033003,      10.7200022, 
	     9.33634567,       8.2426796,      7.40928364,      6.79297686, 
	     6.34373522,      6.01063871,      5.74666357,      5.51208353, 
	     5.27642059,      5.01907587,      4.72886801,      4.40278339, 
	     4.04424238,      3.66118288,      3.26418185,      2.86479568, 
	     2.47422004,       2.1023159,      1.75699866,      1.44394588, 
	     1.16657102,     0.926194072,     0.722342134,     0.553120613, 
	    0.415606558,     0.306226671,     0.221093416,     0.156284735, 
	    0.108059399,    0.0730084404,    0.0481470376,    0.0309553891, 
	   0.0193785131,    0.0117959473,   0.00697184261,    0.0039948686, 
	  0.00221563689,   0.00118740648,  0.000613805256,  0.000305476977, 
	 0.000146079372,  6.69831279e-05,  2.93878602e-05,  1.23086456e-05, 
	 4.90973798e-06,  1.86049226e-06,  6.68009079e-07,  2.26637297e-07, 
	 7.24482589e-08,  2.17554614e-08,  6.11771656e-09, 
/* Z=55 -> Cs */
	              0,    0.0853344053,    0.0934358165,     0.102301657, 
	    0.112003364,     0.122618876,      0.13423337,      0.14693974, 
	    0.160839409,     0.176043004,     0.192671165,     0.210855514, 
	    0.230739504,     0.252479464,     0.276245773,     0.302223951, 
	    0.330616027,     0.361641943,     0.395540982,     0.432573378, 
	    0.473022133,     0.517194748,     0.565425277,     0.618076324, 
	     0.67554152,     0.738247693,     0.806657553,     0.881272256, 
	    0.962634563,       1.0513314,       1.1479975,      1.25331843, 
	     1.36803401,      1.49294245,      1.62890339,      1.77684236, 
	     1.93775415,      2.11270738,      2.30284834,      2.50940514, 
	     2.73369169,       2.9771111,      3.24116015,      3.52743292, 
	     3.83762264,      4.17352581,      4.53704357,      4.93018341, 
	     5.35505867,      5.81389046,      6.30900192,      6.84281921, 
	     7.41786385,      8.03674507,      8.70215416,      9.41684723, 
	     10.1836357,      11.0053635,      11.8848848,      12.8250446, 
	     13.8286409,      14.8983889,      16.0368843,      17.2465572, 
	     18.5296097,      19.8879719,      21.3232193,      22.8365211, 
	     24.4285412,      26.0993671,      27.8484154,      29.6743374, 
	     31.5749245,      33.5469971,      35.5863113,      37.6874695, 
	     39.8437958,       42.047287,      44.2885094,      46.5565605, 
	     48.8390236,      51.1219673,      53.3899612,      55.6261673, 
	     57.8124199,      59.9294281,      61.9569817,      63.8742599, 
	     65.6601868,      67.2938766,      68.7551498,      70.0251083, 
	     71.0867996,      71.9259262,      72.5315857,       72.897049, 
	     73.0205383,      72.9059525,      72.5635376,      72.0104218, 
	     71.2710114,      70.3771515,      69.3680573,      68.2898636, 
	     67.1949234,      66.1406174,      65.1878357,      64.3989792, 
	     63.8356323,      63.5558472,      63.6112175,      64.0437622, 
	     64.8828049,      66.1420517,      67.8169785,      69.8828659, 
	      72.293541,      74.9812393,      77.8575516,      80.8157806, 
	     83.7346725,      86.4835663,      88.9288025,      90.9412308, 
	     92.4043884,      93.2230225,      93.3311539,      92.6992798, 
	     91.3398819,      89.3105545,      86.7141953,      83.6957703, 
	     80.4354095,      77.1378555,      74.0187073,      71.2881622, 
	     69.1334076,      67.7011261,      67.0817413,      67.2971497, 
	     68.2934952,      69.9403992,      72.0372391,      74.3266907, 
	     76.5147705,      78.2959213,      79.3808975,      79.5247879, 
	       78.55233,      76.3774796,      73.0151215,      68.5830231, 
	     63.2935829,      57.4356308,      51.3480186,      45.3874931, 
	     39.8942108,      35.1587372,      31.3941727,      28.7168179, 
	     27.1375313,      26.5649719,      26.8200378,      27.6596375, 
	     28.8064823,       29.980957,      30.9309559,      31.4560776, 
	     31.4236355,      30.7752304,      29.5241337,      27.7448788, 
	     25.5573368,      23.1079693,      20.5508633,      18.0308399, 
	      15.670167,      13.5598831,      11.7558413,      10.2791586, 
	     9.12022209,      8.24520016,      7.60389042,      7.13779354, 
	     6.78750563,      6.49873638,      6.22658205,      5.93793726, 
	     5.61217403,      5.24038601,      4.82359457,      4.37034321, 
	     3.89406562,      3.41057229,      2.93586326,      2.48443198, 
	     2.06810069,      1.69537675,      1.37127006,      1.09747493, 
	    0.872815192,     0.693846822,     0.555525243,     0.451861143, 
	    0.376508713,     0.323249787,     0.286357015,     0.260834724, 
	    0.242548153,     0.228259966,     0.215594977,     0.202956125, 
	    0.189410821,     0.174563527,     0.158426881,     0.141299039, 
	     0.12365289,     0.106040247,    0.0890132934,     0.073064208, 
	   0.0585834607,    0.0458366498,    0.0349580683,    0.0259587727, 
	   0.0187456626,    0.0131475357,   0.00894395076,   0.00589299481, 
	  0.00375498645,   0.00231021619,   0.00137004955, 
/* Z=56 -> Ba */
	              0,    0.0928303152,     0.101608858,      0.11121235, 
	     0.12171752,     0.133208126,     0.145775586,     0.159519613, 
	    0.174548954,     0.190982282,      0.20894888,     0.228589743, 
	    0.250058502,     0.273522466,     0.299163908,      0.32718125, 
	     0.35779044,     0.391226411,     0.427744687,     0.467623025, 
	    0.511163294,     0.558693349,     0.610569179,     0.667176962, 
	    0.728935719,     0.796299398,     0.869760096,     0.949850261, 
	     1.03714621,      1.13227081,      1.23589718,       1.3487519, 
	     1.47161865,      1.60534179,      1.75083065,      1.90906274, 
	     2.08108878,      2.26803589,      2.47111249,      2.69161081, 
	     2.93091369,      3.19049525,      3.47192693,      3.77687979, 
	     4.10712719,      4.46454811,       4.8511281,      5.26896048, 
	     5.72024679,      6.20729446,      6.73251581,      7.29842234, 
	       7.907619,      8.56279659,      9.26672173,      10.0222187, 
	       10.83216,        11.69944,      12.6269522,      13.6175623, 
	     14.6740713,        15.79918,      16.9954433,      18.2652206, 
	     19.6106186,      21.0334263,      22.5350533,      24.1164455, 
	     25.7780037,      27.5195045,      29.3399944,      31.2377014, 
	     33.2099266,      35.2529488,      37.3619194,      39.5307617, 
	     41.7520866,      44.0170937,      46.3155327,      48.6356201, 
	     50.9640541,      53.2859802,       55.585083,      57.8436241, 
	     60.0426369,      62.1620712,      64.1810913,      66.0783997, 
	      67.832634,      69.4228363,      70.8290253,      72.0328217, 
	     73.0181198,      73.7718353,      74.2846985,      74.5520172, 
	     74.5744324,      74.3587036,      73.9182587,      73.2737198, 
	     72.4532242,      71.4924545,       70.434433,      69.3289642, 
	     68.2316895,      67.2027588,      66.3050461,      65.6020279, 
	     65.1551666,      65.0211105,      65.2485504,      65.8750992, 
	     66.9240952,      68.4017944,      70.2949448,      72.5690155, 
	     75.1674576,      78.0119247,      81.0038834,      84.0276108, 
	     86.9545975,      89.6494217,      91.9767761,      93.8095093, 
	     95.0372086,      95.5747757,        95.37043,      94.4123764, 
	     92.7335205,      90.4133911,      87.5768967,      84.3893509, 
	     81.0477066,      77.7682114,      74.7709885,      72.2625351, 
	     70.4174652,      69.3610382,      69.1543274,      69.7837067, 
	     71.1562881,      73.1023407,      75.3852921,      77.7190094, 
	     79.7912064,      81.2911987,      81.9393845,      81.5154343, 
	     79.8822632,      77.0028305,      72.9476166,      67.8915405, 
	     62.0999832,      55.9051476,      49.6747513,      43.7763863, 
	     38.5413094,      34.2318344,      31.0160656,      28.9532108, 
	     27.9911118,      27.9764309,      28.6759911,      29.8063755, 
	     31.0677929,      32.1777115,      32.8999748,      33.0660439, 
	     32.5863419,      31.4511814,      29.7223797,      27.5176086, 
	     24.9904137,      22.3088551,      19.6355629,      17.1113377, 
	     14.8436747,      12.9007902,      11.3109245,      10.0662003, 
	     9.12978077,      8.44501114,      7.94510365,      7.56217098, 
	      7.2346859,      6.91278267,      6.56121063,      6.16005087, 
	     5.70357609,      5.19776344,      4.65703154,      4.10073519, 
	     3.54986501,      3.02427053,      2.54059768,      2.11100459, 
	     1.74262881,      1.43768883,      1.19408071,      1.00629425, 
	    0.866483867,     0.765551448,     0.694122314,     0.643342674, 
	    0.605454981,     0.574147284,      0.54469496,     0.513931632, 
	    0.480094701,     0.442590594,     0.401720911,     0.358403414, 
	    0.313911945,     0.269651771,     0.226979703,     0.187072143, 
	    0.150840476,     0.118890047,    0.0915163681,    0.0687311888, 
	   0.0503098927,    0.0358512774,    0.0248415619,     0.016715169, 
	   0.0109068146,   0.00689128879,   0.00420959434,   0.00248199003, 
	  0.00141001039,  0.000770384911,  0.000404030259, 
/* Z=57 -> La */
	              0,     0.100932971,     0.110439457,     0.120835572, 
	    0.132203683,      0.14463371,     0.158223674,     0.173080504, 
	    0.189320803,     0.207071632,     0.226471484,     0.247671276, 
	     0.27083534,     0.296142638,     0.323787987,     0.353983372, 
	    0.386959404,     0.422966868,     0.462278306,     0.505189955, 
	     0.55202347,     0.603128076,     0.658882737,      0.71969837, 
	    0.786020458,     0.858331621,     0.937154293,      1.02305377, 
	      1.1166414,      1.21857762,      1.32957518,      1.45040345, 
	     1.58189154,      1.72493184,      1.88048482,      2.04958248, 
	     2.23333263,       2.4329226,      2.64962459,      2.88479781, 
	     3.13989472,      3.41646266,      3.71614885,       4.0407033, 
	     4.39198017,      4.77194214,      5.18265963,      5.62631369, 
	     6.10519218,      6.62169027,      7.17830658,      7.77763748, 
	     8.42237091,      9.11527443,      9.85918522,      10.6569948, 
	     11.5116282,      12.4260244,      13.4031048,      14.4457493, 
	     15.5567513,      16.7387829,      17.9943447,      19.3257084, 
	     20.7348652,      22.2234516,      23.7926769,      25.4432468, 
	     27.1752758,      28.9881973,       30.880661,       32.850441, 
	     34.8943214,      37.0080032,      39.1860046,      41.4215469, 
	     43.7064857,      46.0312195,      48.3846436,      50.7540932, 
	     53.1253624,      55.4826889,      57.8088531,      60.0852776, 
	     62.2921944,      64.4088669,      66.4138947,      68.2855835, 
	     70.0023727,      71.5433731,      72.8889313,      74.0213394, 
	     74.9254913,      75.5897293,      76.0065765,      76.1735458, 
	     76.0939331,      75.7774734,      75.2409439,      74.5085831, 
	     73.6123199,      72.5916748,       71.493454,      70.3709717, 
	     69.2829666,      68.2920303,       67.462677,      66.8589172, 
	     66.5415878,      66.5652618,      66.9750977,      67.8035812, 
	     69.0674362,      70.7648621,      72.8733597,      75.3482971, 
	     78.1225586,      81.1073532,      84.1944275,      87.2597122, 
	     90.1684799,      92.7818298,      94.9643936,      96.5928116, 
	     97.5645523,      97.8065567,      97.2829285,      96.0009537, 
	     94.0148315,      91.4263153,      88.3817978,      85.0655594, 
	     81.6891556,      78.4772415,      75.6507034,      73.4081421, 
	     71.9072266,      71.2476883,      71.4578171,      72.4861603, 
	     74.1999054,      76.3908997,      78.7894211,      81.0851288, 
	     82.9536667,      84.0866089,      84.2218628,      83.1713562, 
	     80.8428268,       77.253067,       72.530571,      66.9068451, 
	       60.69664,      54.2688866,      48.0110474,      42.2907677, 
	     37.4190292,      33.6191826,         31.0056,      29.5747204, 
	     29.2095699,      29.6971931,      30.7566395,      32.0736465, 
	      33.337368,      34.2742996,      34.6752014,      34.4121017, 
	      33.444046,      31.8119564,      29.6244659,      27.0375595, 
	     24.2313061,      21.3869381,      18.6669693,      16.2003765, 
	     14.0737619,      12.3287373,      10.9648209,      9.94668674, 
	     9.21417236,      8.69341469,      8.30757809,      7.98594284, 
	     7.67055607,      7.32008505,       6.9109478,      6.43612719, 
	     5.90228939,      5.32591248,      4.72910452,      4.13568497, 
	     3.56792355,      3.04418349,      2.57753348,      2.17526722, 
	     1.83918154,      1.56640184,      1.35053766,      1.18295729, 
	     1.05401254,     0.954086244,     0.874394178,     0.807511926, 
	     0.74764353,     0.690669298,     0.634026706,     0.576482177, 
	    0.517844915,      0.45866558,     0.399951488,     0.342918396, 
	    0.288789958,     0.238648713,     0.193337798,     0.153408483, 
	    0.119107358,    0.0903953463,    0.0669896901,    0.0484208241, 
	   0.0340954289,    0.0233585313,    0.0155486688,    0.0100419549, 
	  0.00628297077,   0.00380226155,   0.00222187932,   0.00125150639, 
	 0.000678222626,  0.000352929899,  0.000175989961, 
/* Z=58 -> Ce */
	              0,     0.109678723,     0.119966492,     0.131212875, 
	    0.143506274,     0.156943038,     0.171628267,     0.187676474, 
	    0.205212533,     0.224372461,     0.245304406,     0.268169761, 
	    0.293144137,     0.320418805,     0.350201726,     0.382719219, 
	    0.418217301,     0.456963301,     0.499247789,     0.545386255, 
	    0.595721066,     0.650623977,     0.710497797,     0.775779486, 
	    0.846942067,     0.924497962,      1.00900149,      1.10105193, 
	     1.20129704,      1.31043613,      1.42922366,      1.55847287, 
	     1.69905984,      1.85192692,      2.01808715,       2.1986289, 
	     2.39471817,      2.60760498,      2.83862615,      3.08921027, 
	      3.3608799,      3.65525794,      3.97406816,      4.31914043, 
	     4.69241095,      5.09592676,      5.53184462,      6.00243044, 
	     6.51006031,      7.05721521,      7.64647961,      8.28053284, 
	     8.96214008,      9.69414425,      10.4794502,      11.3210068, 
	     12.2217894,      13.1847687,      14.2128916,      15.3090353, 
	     16.4759769,      17.7163448,      19.0325661,      20.4268169, 
	     21.9009457,      23.4564114,      25.0942097,      26.8147793, 
	     28.6179218,      30.5026951,      32.4673271,      34.5091095, 
	     36.6242752,      38.8079262,      41.0539055,      43.3547134, 
	     45.7014236,      48.0836143,      50.4893036,      52.9049377, 
	     55.3153954,      57.7040215,      60.0527115,      62.3420601, 
	     64.5515518,      66.6598129,      68.6449585,      70.4849701, 
	     72.1582108,      73.6439667,      74.9230576,      75.9785767, 
	     76.7966232,      77.3670959,       77.684494,      77.7487335, 
	     77.5658722,      77.1488113,      76.5177612,      75.7006607, 
	     74.7332306,      73.6588211,      72.5278854,      71.3970795, 
	     70.3279572,      69.3852158,      68.6345367,      68.1399918, 
	     67.9611511,      68.1499329,      68.7472839,      69.7799377, 
	     71.2574005,      73.1693268,      75.4836197,      78.1453476, 
	     81.0769196,      84.1793976,      87.3353882,      90.4133987, 
	     93.2736435,      95.7751923,      97.7842255,       99.182869, 
	     99.8781891,      99.8107071,      98.9616394,      97.3581543, 
	     95.0759354,      92.2384033,      89.0121078,      85.5981598, 
	     82.2198029,      79.1066437,      76.4765396,      74.5163727, 
	     73.3635025,      73.0896072,      73.6889572,      75.0727005, 
	     77.0704498,      79.4399567,      81.8845139,        84.07724, 
	     85.6901703,      86.4256287,      86.0466156,      84.4029388, 
	     81.4500198,      77.2577591,      72.0080109,      65.9803772, 
	     59.5271301,       53.039753,      46.9104652,      41.4929085, 
	     37.0665665,      33.8092155,      31.7807598,      30.9206505, 
	     31.0590878,      31.9405136,      33.2561569,      34.6811523, 
	     35.9112549,      36.6944122,      36.8535461,      36.2982712, 
	     35.0251999,      33.1079826,      30.6796093,      27.9101257, 
	     24.9831676,      22.0743465,      19.3337364,      16.8740406, 
	     14.7647972,      13.0323839,      11.6648521,      10.6201773, 
	     9.83639431,      9.24208355,      8.76591492,       8.3443222, 
	     7.92676401,      7.47844076,      6.98067713,      6.42941236, 
	     5.83238173,      5.20559454,      4.56966448,      3.94643211, 
	     3.35618496,       2.8156445,      2.33675981,      1.92625558, 
	      1.5858134,      1.31272805,      1.10086644,     0.941763461, 
	    0.825715899,      0.74276334,     0.683486164,     0.639583051, 
	    0.604224622,     0.572202146,      0.53990835,     0.505191922, 
	    0.467131615,     0.425768346,     0.381828636,     0.336463332, 
	     0.29101783,     0.246843785,     0.205156475,     0.166937575, 
	    0.132881284,     0.103378028,    0.0785297081,    0.0581884049, 
	    0.042010501,    0.0295179337,    0.0201591961,    0.0133639984, 
	  0.00858728494,   0.00534041505,   0.00320921722,   0.00186034909, 
	  0.00103845296,  0.000557134685,  0.000286717783, 
/* Z=59 -> Pr */
	              0,     0.119143568,     0.130272105,     0.142433077, 
	    0.155721217,     0.170239791,     0.186101347,     0.203428507, 
	    0.222354844,     0.243025884,     0.265599996,     0.290249646, 
	    0.317162484,     0.346542567,     0.378611922,     0.413611889, 
	    0.451804698,     0.493475258,     0.538932979,     0.588513553, 
	    0.642581344,      0.70153147,     0.765792012,     0.835826993, 
	    0.912138462,     0.995269954,      1.08580899,      1.18439054, 
	     1.29170036,      1.40847814,      1.53552163,      1.67368984, 
	     1.82390749,      1.98716891,      2.16454172,      2.35717154, 
	     2.56628585,      2.79319835,      3.03931332,       3.3061285, 
	     3.59524107,       3.9083488,      4.24725485,      4.61386967, 
	     5.01021338,      5.43841791,      5.90072536,      6.39948988, 
	     6.93717337,      7.51634312,      8.13966751,      8.80990505, 
	      9.5298996,      10.3025627,      11.1308641,      12.0178041, 
	     12.9663992,      13.9796495,      15.0605097,      16.2118511, 
	     17.4364204,      18.7367859,      20.1152878,      21.5739822, 
	      23.114563,      24.7382927,      26.4459209,      28.2375984, 
	     30.1127834,      32.0701408,      34.1074371,      36.2214508, 
	     38.4078407,      40.6610718,      42.9742889,       45.339241, 
	     47.7461967,      50.1838837,      52.6394424,      55.0984116, 
	     57.5447655,      59.9609413,      62.3279877,      64.6257019, 
	     66.8328629,      68.9275208,      70.8873825,      72.6902084, 
	     74.3143997,      75.7395401,      76.9470825,      77.9210968, 
	     78.6490326,      79.1225586,      79.3383255,      79.2988358, 
	     79.0131073,      78.4973221,      77.7752914,      76.8787079, 
	     75.8471527,      74.7277527,      73.5745239,      72.4472885, 
	     71.4101562,      70.5295715,       69.871933,      69.5008621, 
	     69.4740829,      69.8401794,       70.635231,        71.87957, 
	     73.5748444,      75.7016296,      78.2178116,      81.0579681, 
	     84.1339798,      87.3370438,      90.5411682,      93.6082611, 
	     96.3946304,      98.7587357,      100.569908,       101.71756, 
	     102.120163,      101.733528,       100.55751,      98.6402359, 
	     96.0794678,      93.0201874,       89.648262,      86.1799393, 
	      82.847702,      79.8829269,      77.4967117,      75.8602524, 
	     75.0866394,      75.2159729,      76.2057648,      77.9281464, 
	     80.1750031,      82.6713867,      85.0966415,      87.1118088, 
	     88.3909836,      88.6535721,      87.6941528,      85.4064102, 
	     81.7982712,      76.9959106,      71.2356567,      64.8439331, 
	     58.2069626,      51.7331848,      45.8123894,      40.7761612, 
	      36.864399,      34.2019806,      32.7885666,      32.5029335, 
	     33.1211777,      34.3462334,      35.8447113,      37.2859383, 
	     38.3780441,      38.8966484,      38.7028961,      37.7496338, 
	     36.0760193,      33.7926865,      31.0604286,      28.0659885, 
	     24.9982986,      22.0280228,      19.2923145,      16.8858223, 
	     14.8579226,      13.2155094,      11.9300251,      10.9471903, 
	      10.197793,      9.60808563,      9.10865021,      8.64099026, 
	     8.16154003,      7.64315081,      7.07442236,      6.45741892, 
	     5.80438375,      5.13405132,      4.46805191,      3.82779574, 
	      3.2320478,      2.69531107,      2.22699857,      1.83130097, 
	     1.50761271,      1.25134754,      1.05497491,     0.909125924, 
	    0.803642869,     0.728483081,     0.674422204,     0.633535206, 
	    0.599462211,     0.567486048,     0.534459531,     0.498626709, 
	     0.45937866,     0.416981578,     0.372304857,     0.326570183, 
	    0.281135887,     0.237323165,     0.196287602,     0.158935085, 
	    0.125878707,    0.0974314436,     0.073627606,    0.0542657562, 
	   0.0389646105,    0.0272242147,    0.0184853226,    0.0121813631, 
	  0.00777922291,   0.00480712252,    0.0028697364,   0.00165221572, 
	 0.000915754063,   0.00048770319,  0.000249075121, 
/* Z=60 -> Nd */
	              0,     0.129378513,     0.141410962,      0.15455471, 
	    0.168911234,     0.184591115,     0.201714784,     0.220413461, 
	    0.240829974,     0.263119817,     0.287452281,     0.314011425, 
	    0.342997521,     0.374628365,     0.409140646,      0.44679153, 
	    0.487860352,     0.532650471,     0.581491053,     0.634739161, 
	    0.692782044,     0.756039441,     0.824965954,      0.90005368, 
	    0.981835246,      1.07088649,      1.16782975,      1.27333713, 
	     1.38813365,      1.51300144,      1.64878285,      1.79638457, 
	     1.95678186,      2.13102245,      2.32023025,      2.52561069, 
	     2.74845386,      2.99013901,      3.25213933,      3.53602481, 
	     3.84346676,      4.17624092,      4.53623009,      4.92542696, 
	     5.34593678,      5.79997587,       6.2898736,      6.81806993, 
	     7.38711309,       7.9996562,      8.65844727,      9.36632633, 
	     10.1262074,      10.9410696,        11.81394,      12.7478647, 
	     13.7458973,      14.8110542,       15.946291,      17.1544571, 
	     18.4382553,      19.8001823,      21.2424755,      22.7670498, 
	     24.3754215,      26.0686359,      27.8471756,      29.7108746, 
	     31.6588211,      33.6892548,      35.7994576,      37.9856606, 
	     40.2429199,       42.565033,      44.9444199,      47.3720512, 
	     49.8373756,      52.3282394,      54.8308983,      57.3299828, 
	     59.8085556,      62.2481766,      64.6290588,      66.9302368, 
	     69.1298523,      71.2054443,      73.1343765,      74.8943253, 
	     76.4638214,      77.8228989,      78.9537888,      79.8417053, 
	     80.4756317,      80.8491669,       80.961319,       80.817337, 
	     80.4293442,      79.8169708,      79.0077057,      78.0370483, 
	     76.9484406,      75.7927322,      74.6273499,      73.5150681, 
	     72.5222244,      71.7166138,      71.1649017,      70.9296646, 
	     71.0661926,      71.6190872,      72.6189041,      74.0789185, 
	     75.9923859,      78.3303375,      81.0403671,      84.0464325, 
	      87.250061,      90.5329895,      93.7613678,      96.7914963, 
	     99.4770584,      101.677422,      103.266792,       104.14357, 
	     104.239372,       103.52684,       102.02562,       99.805542, 
	     96.9865189,      93.7344589,       90.253067,      86.7715378, 
	     83.5286026,      80.7539062,      78.6479721,      77.3624802, 
	      76.982811,      77.5148392,      78.8778076,      80.9047852, 
	     83.3515472,       85.913765,      88.2517014,      90.0204239, 
	     90.9029007,      90.6425171,      89.0716171,      86.1325836, 
	     81.8886032,      76.5223541,      70.3220215,      63.6555939, 
	     56.9356194,      50.5782776,      44.9610634,      40.3841705, 
	      37.040184,       34.996048,      34.1894836,      34.4405098, 
	     35.4763641,      36.9663506,      38.5618362,      39.9359207, 
	      40.817627,      41.0165977,      40.4358215,       39.071949, 
	     37.0045547,      34.3769455,      31.3721581,      28.1876717, 
	     25.0122051,      22.0070992,      19.2937317,      16.9475231, 
	     14.9981031,      13.4345541,      12.2142401,      11.2735205, 
	     10.5387478,      9.93622684,      9.40014648,      8.87800503, 
	     8.33338642,      7.74634457,      7.11186171,      6.43699217, 
	     5.73729849,      5.03314638,      4.34629393,      3.69707942, 
	     3.10237026,      2.57431555,      2.11984968,      1.74082518, 
	     1.43462634,      1.19509184,      1.01359022,     0.880104661, 
	    0.784225166,     0.715968847,     0.666391969,     0.627983153, 
	    0.594853759,     0.562754989,     0.528963923,     0.492078841, 
	    0.451764375,     0.408478409,     0.363206953,     0.317224354, 
	    0.271890849,      0.22849232,      0.18812488,      0.15162231, 
	    0.119522683,    0.0920686424,    0.0692344233,     0.050772015, 
	   0.0362683944,    0.0252063051,    0.0170218702,    0.0111539168, 
	  0.00708170794,    0.0043498138,   0.00258058822,   0.00147617003, 
	 0.000812712766,  0.000429825159,  0.000217935289, 
/* Z=61 -> Pm */
	              0,     0.140446872,     0.153451085,     0.167650819, 
	    0.183154732,      0.20008114,     0.218558863,     0.238728195, 
	     0.26074174,     0.284765601,     0.310980469,     0.339582801, 
	     0.37078616,     0.404822737,     0.441944659,     0.482425869, 
	    0.526563644,     0.574680686,     0.627126992,     0.684282005, 
	    0.746556997,     0.814397454,     0.888285577,     0.968743145, 
	      1.0563345,      1.15166938,      1.25540614,      1.36825573, 
	     1.49098408,      1.62441719,      1.76944387,      1.92702031, 
	     2.09817386,      2.28400731,      2.48570299,      2.70452785, 
	     2.94183612,      3.19907475,      3.47778726,      3.77961731, 
	     4.10631227,      4.45972681,      4.84182549,      5.25468397, 
	     5.70049143,      6.18155098,      6.70027828,      7.25919962, 
	     7.86094856,      8.50826073,      9.20396519,      9.95097542, 
	     10.7522764,      11.6109066,      12.5299473,      13.5124846, 
	      14.561595,      15.6803074,      16.8715649,      18.1381836, 
	     19.4828072,      20.9078445,      22.4154129,      24.0072651, 
	     25.6847248,      27.4485874,      29.2990494,      31.2355995, 
	     33.2569275,      35.3608208,      37.5440483,      39.8022614, 
	     42.1298866,      44.5200195,      46.9643326,      49.4529953, 
	     51.9746017,      54.5161285,      57.0629196,      59.5986977, 
	     62.1056213,      64.5643921,       66.954422,      69.2540283, 
	     71.4407654,      73.4917297,      75.3840714,      77.0954514, 
	     78.6046829,      79.8923798,      80.9417114,      81.7392044, 
	     82.2755508,       82.546463,      82.5534897,      82.3047791, 
	     81.8157349,        81.10952,      80.2173843,      79.1787109, 
	     78.0407257,      76.8579254,      75.6909866,      74.6053467, 
	     73.6692581,      72.9514389,      72.5182877,      72.4307861, 
	     72.7411499,       73.489357,      74.6997833,      76.3780594, 
	      78.508522,       81.052269,      83.9463806,      87.1042099, 
	     90.4171753,      93.7580948,      96.9860535,      99.9529037, 
	     102.511047,      104.522369,      105.867622,      106.455971, 
	     106.233803,      105.191994,      103.370972,       100.86274, 
	     97.8091507,      94.3960876,      90.8433533,      87.3905487, 
	     84.2795029,      81.7344818,      79.9416122,      79.0294037, 
	     79.0524063,      79.9800644,      81.6924515,      83.9842377, 
	     86.5773773,      89.1421509,       91.325119,      92.7817993, 
	     93.2108383,      92.3861923,      90.1835785,      86.5979538, 
	     81.7493591,      75.8759079,      69.3138504,      62.4663315, 
	     55.7638855,      49.6208878,      44.3929329,      40.3402596, 
	     37.6017532,      36.1830063,      35.9598923,      36.6972427, 
	      38.079998,      39.7524376,      41.3601379,      42.5889626, 
	     43.1962128,      43.0304832,      42.0386162,      40.2602959, 
	     37.8123016,      34.8657951,      31.6203995,      28.2787762, 
	     25.0246754,      22.0066662,      19.3284473,      17.0458736, 
	     15.1698189,      13.6735411,      12.5028467,      11.5873642, 
	     10.8513756,      10.2230644,      9.64138031,      9.06027508, 
	     8.45032215,      7.79814577,      7.10418653,      6.37944746, 
	     5.64179802,      4.91236782,      4.21237707,      3.56066346, 
	     2.97198844,      2.45612788,      2.01766062,      1.65631831, 
	      1.3677386,      1.14445627,     0.976986051,     0.854870677, 
	    0.767603695,     0.705366671,     0.659555197,      0.62309283, 
	     0.59055531,     0.558138549,     0.523510516,      0.48558861, 
	    0.444277614,     0.400198758,     0.354432523,     0.308290273, 
	    0.263124287,     0.220180437,     0.180494636,     0.144830823, 
	    0.113656253,    0.0871485993,    0.0652276427,    0.0476041064, 
	   0.0338378474,    0.0233978648,    0.0157180876,    0.0102440687, 
	  0.00646781549,   0.00394985452,   0.00232932949,   0.00132420578, 
	 0.000724373269,  0.000380554528,  0.000191619431, 
/* Z=62 -> Sm */
	              0,     0.152417436,     0.166466266,     0.181800589, 
	    0.198536649,     0.216800973,     0.236731291,      0.25847742, 
	    0.282202363,     0.308083415,     0.336313337,      0.36710161, 
	    0.400675893,     0.437283427,      0.47719276,     0.520695329, 
	    0.568107426,     0.619771898,     0.676060736,     0.737376809, 
	    0.804156542,     0.876872361,     0.956035495,      1.04219866, 
	     1.13595915,      1.23796201,      1.34890354,      1.46953452, 
	      1.6006639,      1.74316299,      1.89796889,      2.06608868, 
	     2.24860382,      2.44667411,      2.66154242,      2.89453769, 
	     3.14708138,      3.42068911,      3.71697688,      4.03766298, 
	     4.38457298,      4.75964165,      5.16491461,      5.60255337, 
	     6.07483101,      6.58413792,       7.1329751,      7.72395468, 
	     8.35979366,      9.04331017,       9.7774086,      10.5650749, 
	     11.4093599,      12.3133602,      13.2801971,       14.312993, 
	     15.4148407,      16.5887661,      17.8376923,      19.1643906, 
	     20.5714264,      22.0611038,      23.6354008,      25.2958927, 
	     27.0436783,      28.8792934,      30.8026104,      32.8127518, 
	     34.9079819,      37.0856018,      39.3418388,      41.6717377, 
	     44.0690613,      46.5261803,      49.0339966,      51.5818481, 
	     54.1574707,      56.7469368,      59.3346863,      61.9035339, 
	     64.4347534,      66.9082031,      69.3025284,      71.5953903, 
	     73.7638092,      75.7845459,      77.6345978,      79.2917557, 
	     80.7352295,      81.9463806,      82.9094696,      83.6124954, 
	     84.0480576,      84.2141495,      84.1150284,      83.7619019, 
	      83.173584,      82.3768921,      81.4069061,      80.3068542, 
	     79.1277542,      77.9275818,      76.7700653,      75.7230377, 
	     74.8562622,      74.2389297,      73.9366531,      74.0082245, 
	     74.5021362,      75.4531174,      76.8787231,      78.7763901, 
	     81.1210403,      83.8635406,      86.9303055,      90.2241898, 
	     93.6268692,       97.002861,      100.205162,      103.082352, 
	     105.487038,       107.28521,      108.365921,      108.650803, 
	     108.102539,      106.731483,      104.599663,      101.821442, 
	     98.5601349,      95.0203094,      91.4358521,      88.0539627, 
	     85.1162262,      82.8378296,      81.3867264,      80.8648148, 
	     81.2930832,      82.6028748,      84.6347961,      87.1464615, 
	       89.82901,       92.331871,      94.2937241,      95.3771133, 
	     95.3031235,      93.8824081,      91.0389252,      86.8230515, 
	     81.4121323,      75.0973969,      68.2581635,      61.3254547, 
	     54.7388649,      48.9013939,      44.1376076,      40.6602402, 
	     38.5495071,      37.7479286,       38.071312,      39.2343292, 
	     40.8871994,      42.6582031,      44.1963577,      45.2084694, 
	     45.4862518,      44.9207039,      43.5031548,      41.3143234, 
	     38.5042305,      35.2667122,      31.8124981,      28.3443794, 
	     25.0371418,       22.023983,      19.3898411,      17.1713181, 
	     15.3619499,       13.921298,      12.7860374,      11.8813963, 
	     11.1315374,      10.4678907,      9.83490658,       9.1931057, 
	     8.51967716,      7.80709505,      7.06035709,      6.29347849, 
	     5.52578878,       4.7784934,       4.0717907,      3.42273712, 
	      2.8438921,      2.34271169,      1.92157984,      1.57832909, 
	     1.30709553,      1.09934795,     0.944956422,      0.83319056, 
	    0.753569186,      0.69651562,     0.653803527,     0.618800104, 
	    0.586533248,     0.553618073,       0.5180825,     0.479130745, 
	    0.436878413,     0.392085463,     0.345907807,     0.299680263, 
	    0.254738599,     0.212284625,     0.173293263,     0.138460413, 
	    0.108185999,    0.0825871229,    0.0615337789,    0.0446999185, 
	   0.0316220745,    0.0217584893,    0.0145429103,   0.00942870695, 
	  0.00592091726,    0.0035956935,   0.00210821885,   0.00119132712, 
	 0.000647635607,  0.000338044687,  0.000169073595, 
/* Z=63 -> Eu */
	              0,      0.16536583,     0.180537641,     0.197091043, 
	    0.215150326,     0.234850675,     0.256339252,     0.279775977, 
	    0.305334866,     0.333205014,     0.363591909,      0.39671877, 
	    0.432828128,     0.472183138,     0.515069544,     0.561797261, 
	     0.61270231,     0.668149173,     0.728532493,     0.794279814, 
	    0.865853965,     0.943755567,      1.02852607,      1.12075055, 
	     1.22106063,      1.33013844,      1.44871938,      1.57759595, 
	     1.71762168,      1.86971486,      2.03486228,       2.2141242, 
	       2.408638,      2.61962199,      2.84838057,      3.09630847, 
	     3.36489415,      3.65572476,      3.97048974,      4.31098413, 
	     4.67911196,      5.07688951,      5.50644588,      5.97002649, 
	     6.46999168,      7.00881577,      7.58908701,      8.21350193, 
	     8.88485909,      9.60605526,      10.3800697,      11.2099581, 
	     12.0988283,      13.0498276,      14.0661192,      15.1508455, 
	     16.3071098,      17.5379238,      18.8461761,      20.2345772, 
	     21.7056007,      23.2614326,      24.9038868,      26.6343422, 
	     28.4536514,      30.3620567,      32.3590927,      34.4434776, 
	     36.6130295,      38.8645287,       41.193634,      43.5947533, 
	     46.0609627,      48.5838737,      51.1535912,      53.7586136, 
	     56.3857841,      59.0202827,      61.6456299,      64.2437363, 
	     66.7950211,      69.2785187,      71.6721573,      73.9529877, 
	     76.0975723,      78.0824203,      79.8845062,      81.4818497, 
	     82.8542023,      83.9838257,      84.8562317,      85.4610901, 
	     85.7930527,      85.8525925,      85.6468201,      85.1901855, 
	     84.5050049,      83.6218643,      82.5796814,        81.42556, 
	     80.2141876,      79.0068817,      77.8701935,      76.8740158, 
	     76.0892334,       75.585022,      75.4256516,      75.6671143, 
	      76.353569,      77.5138092,      79.1580048,      81.2748337, 
	      83.829422,      86.7621231,      89.9886703,      93.4015732, 
	     96.8732376,      100.260674,      103.411751,      106.173126, 
	     108.399117,      109.961502,      110.759331,      110.728363, 
	     109.848991,      108.152161,      105.722076,      102.695442, 
	     99.2562103,      95.6261368,      92.0508194,      88.7821121, 
	     86.0578308,      84.0803375,      82.9957733,      82.8761978, 
	     83.7066498,      85.3791733,      87.6952286,      90.3772659, 
	     93.0891876,      95.4645462,      97.1402588,       97.792511, 
	     97.1713028,      95.1296921,      91.6439362,       86.821846, 
	     80.8974304,      74.2117691,      67.1814041,      60.2573967, 
	     53.8792877,      48.4293976,      44.1928177,      41.3283272, 
	     39.8537903,      39.6483498,      40.4707489,      41.9914703, 
	     43.8340225,      45.6196709,        47.00951,      47.7385826, 
	     47.6380539,      46.6437225,      44.7911758,      42.1996841, 
	     39.0484161,      35.5489807,      31.9182301,      28.3547077, 
	     25.0208855,      22.0325241,      19.4550323,       17.306097, 
	     15.5631132,      14.1736364,      13.0670652,      12.1659622, 
	     11.3957653,      10.6921043,      10.0054045,      9.30281544, 
	     8.56787014,      7.79841089,      7.00341129,      6.19932127, 
	     5.40642166,       4.6456008,       3.9357748,       3.2920835, 
	     2.72486067,      2.23930407,       1.8357197,      1.51019144, 
	     1.25551331,      1.06223989,     0.919730365,     0.817089617, 
	    0.743942976,     0.691010892,     0.650477588,     0.616167665, 
	    0.583560884,     0.549681306,     0.512900233,      0.47268945, 
	    0.429354697,     0.383774847,     0.337162971,     0.290861249, 
	    0.246175826,     0.204253614,     0.166000694,     0.132039353, 
	    0.102698877,     0.078034237,    0.0578655154,    0.0418307222, 
	   0.0294445027,    0.0201560203,    0.0134004913,   0.00864052866, 
	  0.00539530301,   0.00325733447,   0.00189826428,   0.00106594805, 
	 0.000575700367,  0.000298463507,  0.000148227424, 
/* Z=64 -> Gd */
	              0,     0.179383636,     0.195763588,     0.213627815, 
	    0.233109131,     0.254351944,     0.277513355,     0.302764058, 
	    0.330289602,      0.36029163,     0.392989039,     0.428619623, 
	    0.467441499,     0.509734631,     0.555802822,     0.605975509, 
	    0.660609603,     0.720091879,      0.78484112,     0.855310678, 
	    0.931990921,      1.01541209,      1.10614717,      1.20481515, 
	     1.31208348,      1.42867279,      1.55535924,      1.69297886, 
	     1.84243119,      2.00468349,      2.18077421,      2.37181854, 
	     2.57901096,      2.80363107,      3.04704642,      3.31071877, 
	     3.59620619,      3.90516853,      4.23937035,      4.60068464, 
	     4.99109602,      5.41270304,       5.8677187,      6.35847378, 
	     6.88741446,      7.45710087,      8.07020664,      8.72951126, 
	     9.43789387,      10.1983271,      11.0138626,      11.8876181, 
	     12.8227577,      13.8224726,      14.8899565,       16.028368, 
	     17.2408047,      18.5302582,      19.8995686,      21.3513699, 
	     22.8880272,      24.5115871,      26.2236843,      28.0254726, 
	     29.9175415,      31.8998127,      33.9714508,      36.1307564, 
	     38.3750534,      40.7005844,      43.1024055,      45.5742607, 
	     48.1084976,      50.6959686,      53.3259468,      55.9860573, 
	     58.6622581,      61.3388214,      63.9983559,      66.6218796, 
	     69.1889572,      71.6778564,      74.0658188,      76.3293686, 
	     78.4446869,      80.3881149,      82.1367188,      83.6688766, 
	     84.9650497,      86.0085297,      86.7862778,      87.2897949, 
	     87.5159683,      87.4679413,      87.1558456,       86.597496, 
	     85.8188629,      84.8543472,      83.7467804,      82.5470657, 
	     81.3134918,      80.1105499,      79.0073395,      78.0754623, 
	      77.386528,      77.0091629,      77.0057297,      77.4288025, 
	     78.3175354,       79.694191,      81.5609207,      83.8971481, 
	     86.6577759,      89.7724609,      93.1462097,      96.6614685, 
	     100.181892,      103.557739,      106.632957,      109.253593, 
	     111.277237,      112.583191,      113.082283,      112.726044, 
	     111.514061,      109.498901,      106.787521,      103.538811, 
	     99.9565887,      96.2780762,      92.7581787,      89.6502457, 
	        87.1847,      85.5471649,      84.8581161,      85.1563492, 
	     86.3883286,      88.4053802,      90.9697723,      93.7703857, 
	     96.4470978,      98.6223145,      99.9369736,       100.08741, 
	     98.8591156,      96.1534653,      92.0037231,       86.577919, 
	     80.1673126,      73.1608658,      66.0079651,      59.1731682, 
	      53.087925,      48.1050453,      44.4615288,      42.2546349, 
	     41.4342957,      41.8131332,      43.0924988,      44.9010277, 
	     46.8401718,      48.5305061,      49.6525269,      49.9770012, 
	     49.3816719,      47.8534126,        45.47715,      42.4145126, 
	     38.8762741,      35.0929108,      31.2871895,      27.6519203, 
	     24.3346252,      21.4299583,      18.9793472,      16.9766674, 
	     15.3782101,      14.1149321,      13.1051426,      12.2660294, 
	     11.5229006,      10.8155451,      10.1015434,      9.35682392, 
	     8.57401276,      7.75924063,      6.92815542,      6.10176134, 
	     5.30260086,      4.55162764,      3.86595178,      3.25750351, 
	     2.73254919,      2.29191113,      1.93171132,      1.64444089, 
	     1.42016995,      1.24774134,      1.11583078,      1.01379788, 
	    0.932291806,     0.863611758,     0.801846981,     0.742837608, 
	    0.684005439,     0.624101877,     0.562916994,     0.500983417, 
	    0.439299554,     0.379088908,     0.321604013,     0.267977923, 
	    0.219122663,     0.175670937,     0.137955576,     0.106020309, 
	   0.0796543583,    0.0584432371,    0.0418280885,    0.0291665737, 
	   0.0197892711,    0.0130470507,   0.00834662933,   0.00517332554, 
	  0.00310167833,   0.00179582275,   0.00100231241,  0.000538280583, 
	 0.000277607614,  0.000137207666,  6.48501591e-05, 
/* Z=65 -> Tb */
	              0,      0.19452709,     0.212203711,     0.231474027, 
	    0.252479881,     0.275375426,     0.300328255,       0.3275204, 
	    0.357149661,     0.389430761,     0.424596876,     0.462900996, 
	    0.504617631,     0.550044537,     0.599504411,     0.653346956, 
	    0.711951196,     0.775727153,     0.845118999,     0.920606792, 
	     1.00270987,      1.09198928,      1.18905091,      1.29454875, 
	     1.40918791,       1.5337286,      1.66898954,      1.81585169, 
	     1.97526228,      2.14823914,      2.33587432,      2.53933859, 
	     2.75988626,      2.99885893,      3.25768995,      3.53790855, 
	     3.84114456,      4.16913128,      4.52371025,      4.90683317, 
	     5.32056618,      5.76709032,        6.248703,      6.76782036, 
	     7.32697392,      7.92880917,      8.57608318,      9.27165604, 
	      10.018486,      10.8196173,      11.6781693,      12.5973158, 
	     13.5802736,      14.6302681,       15.750514,      16.9441795, 
	     18.2143459,      19.5639668,      20.9958172,      22.5124435, 
	     24.1160908,      25.8086338,      27.5915127,      29.4656372, 
	     31.4313011,      33.4880867,      35.6347618,      37.8691673, 
	      40.188118,      42.5872917,         45.0611,      47.6026039, 
	     50.2034073,      52.8535576,      55.5414848,      58.2539406, 
	      60.975975,      63.6909447,      66.3805618,      69.0249634, 
	     71.6028976,      74.0918961,      76.4685822,       78.708992, 
	     80.7890396,      82.6849976,       84.374115,      85.8352966, 
	     87.0498199,      88.0021896,      88.6809616,      89.0796356, 
	     89.1975098,      89.0405197,      88.6219788,      87.9631805, 
	     87.0938034,      86.0521088,      84.8847885,      83.6464844, 
	     82.3989105,      81.2094803,      80.1495209,      79.2919769, 
	     78.7086639,      78.4671478,      78.6272888,      79.2376251, 
	     80.3317184,      81.9247131,      84.0102615,      86.5581512, 
	     89.5128937,      92.7934113,      96.2942505,      99.8882828, 
	      103.43119,       106.76757,      109.738701,       112.19149, 
	      113.98848,      115.018044,      115.204308,      114.515823, 
	     112.972237,      110.648041,      107.672623,       104.22599, 
	     100.529953,      96.8346863,       93.401268,      90.4810791, 
	     88.2935944,      87.0044785,      86.7060547,      87.4025116, 
	     89.0018921,      91.3166428,      94.0735855,      96.9334335, 
	     99.5188065,      101.448517,      102.375122,      102.021797, 
	     100.214478,      96.9051208,      92.1830139,      86.2718506, 
	     79.5120773,      72.3295898,        65.19384,      58.5696754, 
	     52.8685303,      48.4048042,      45.3629951,      43.7799606, 
	     43.5446167,      44.4152451,      46.0518799,      48.0591812, 
	     50.0338173,      51.6098061,      52.4959526,      52.5010796, 
	     51.5447273,      49.6535606,      46.9454842,      43.6051102, 
	     39.8549232,      35.9262924,      32.0339737,      28.3567066, 
	     25.0250778,      22.1169224,      19.6593933,      17.6362038, 
	     15.9982376,      14.6756001,      13.5894527,       12.662303, 
	     11.8258944,      11.0263224,      10.2264156,      9.40573215, 
	     8.55872059,       7.6916914,      6.81919622,      5.96037006, 
	     5.13561535,      4.36390352,      3.66082239,      3.03739619, 
	     2.49960995,      2.04852772,      1.68085885,      1.38981819, 
	     1.16613626,     0.999093235,     0.877476335,      0.79039216, 
	    0.727892935,     0.681402743,      0.64395237,     0.610245824, 
	    0.576592624,     0.540742278,     0.501656532,     0.459252387, 
	    0.414140016,     0.367376179,     0.320245385,     0.274077624, 
	    0.230105981,     0.189365134,     0.152629048,     0.120384119, 
	   0.0928328484,     0.069921717,    0.0513862744,    0.0368062332, 
	   0.0256635267,    0.0173972063,    0.0114504881,   0.00730680767, 
	   0.0045136963,   0.00269489922,   0.00155247445,  0.000861399167, 
	 0.000459480711,  0.000235155225,  0.000115229261, 
/* Z=66 -> Dy */
	              0,     0.210925415,     0.229997009,     0.250779212, 
	    0.273423374,     0.298094004,     0.324969828,     0.354244858, 
	    0.386129797,     0.420853257,     0.458663404,     0.499829322, 
	    0.544642806,     0.593420208,     0.646504164,     0.704265893, 
	    0.767107308,     0.835463107,     0.909803748,      0.99063766, 
	     1.07851422,      1.17402673,      1.27781534,      1.39057076, 
	      1.5130372,      1.64601636,      1.79037082,       1.9470284, 
	      2.1169858,      2.30131292,      2.50115728,      2.71774769, 
	     2.95239949,      3.20651865,      3.48160529,      3.77925873, 
	     4.10118151,      4.44918251,      4.82518101,      5.23120975, 
	     5.66941595,      6.14206457,      6.65153933,      7.20033979, 
	     7.79108381,      8.42650127,      9.10942936,      9.84281063, 
	     10.6296768,      11.4731417,      12.3763876,      13.3426418, 
	     14.3751621,      15.4772043,      16.6519966,      17.9027004, 
	     19.2323723,      20.6439171,      22.1400337,      23.7231579, 
	     25.3953876,      27.1584244,      29.0134811,      30.9612007, 
	     33.0015564,      35.1337662,      37.3561707,      39.6661263, 
	     42.0599136,      44.5325966,      47.0779381,      49.6882782, 
	     52.3544426,      55.0656586,      57.8094978,      60.5718155, 
	     63.3367577,      66.0867615,      68.8026428,      71.4636917, 
	     74.0478592,      76.5320053,      78.8921661,      81.1040039, 
	      83.143219,      84.9861374,      86.6103439,       87.995369, 
	     89.1234894,      89.9805908,      90.5569839,      90.8483582, 
	     90.8565979,      90.5906067,       90.067009,      89.3106995, 
	     88.3551636,      87.2425613,      86.0234299,      84.7560654, 
	      83.505455,      82.3417053,      81.3380508,      80.5683594, 
	     80.1041946,      80.0114899,      80.3470306,      81.1546707, 
	     82.4617691,      84.2757187,      86.5811691,      89.3378677, 
	     92.4796448,      95.9146652,      99.5271835,      103.180946, 
	     106.724388,       109.99749,       112.84008,      115.101501, 
	     116.650818,       117.38726,      117.249847,      116.225639, 
	     114.355461,      111.736389,      108.520264,      104.907547, 
	     101.136642,       97.468483,      94.1674347,      91.4795151, 
	     89.6095581,      88.6995926,      88.8104401,      89.9090347, 
	     91.8634567,      94.4470291,      97.3523636,      100.214828, 
	     102.643951,      104.260239,      104.733833,      103.820778, 
	     101.392906,      97.4571381,      92.1614304,      85.7857361, 
	     78.7180099,      71.4172287,      64.3671722,      58.0259705, 
	     52.7774086,      48.8900986,      46.4896622,      45.5478935, 
	     45.8901253,      47.2200241,      49.1581078,      51.2886238, 
	     53.2081604,      54.5694809,      55.1151009,      54.6969566, 
	     53.2809906,      50.9376373,      47.8211975,      44.1421623, 
	     40.1370316,      36.0397224,      32.0578423,      28.3559017, 
	     25.0462437,      22.1874065,      19.7886791,      17.8191833, 
	     16.2195187,      14.9141016,      13.8226585,      12.8697147, 
	     11.9914417,       11.139678,      10.2832994,       9.4074173, 
	     8.51098824,       7.6034894,      6.70123148,      5.82381058, 
	     4.99101639,      4.22042322,      3.52574205,      2.91591549, 
	     2.39487791,      1.96184647,       1.6120038,      1.33741331, 
	     1.12803841,     0.972747147,     0.860216677,     0.779677927, 
	    0.721470654,     0.677403212,     0.640930057,     0.607172847, 
	    0.572820961,     0.535945594,     0.495763063,     0.452375859, 
	    0.406515449,     0.359303296,     0.312042534,     0.266046017, 
	    0.222504839,     0.182396367,     0.146430552,     0.115030169, 
	   0.0883399695,    0.0662583262,    0.0484844856,    0.0345743224, 
	   0.0239976868,    0.0161916986,    0.0106054926,    0.0067337607, 
	  0.00413818331,   0.00245746016,   0.00140782644,  0.000776638044, 
	 0.000411788846,  0.000209436752,  0.000101963349, 
/* Z=67 -> Ho */
	              0,     0.228673086,     0.249244317,      0.27165091, 
	    0.296054363,     0.322630078,      0.35156855,     0.383076429, 
	    0.417378038,     0.454716742,     0.495356351,     0.539582968, 
	    0.587706625,     0.640063107,      0.69701612,     0.758959413, 
	      0.8263188,     0.899555027,     0.979165971,      1.06568968, 
	     1.15970695,      1.26184475,      1.37277913,      1.49323893, 
	     1.62400866,      1.76593351,      1.91992164,      2.08694911, 
	     2.26806402,      2.46439028,      2.67713213,      2.90757871, 
	     3.15710831,      3.42719173,      3.71939874,      4.03540039, 
	     4.37697268,      4.74600172,      5.14448643,      5.57454014, 
	     6.03839302,       6.5383954,      7.07701445,      7.65683603, 
	     8.28056145,      8.95100403,      9.67108345,      10.4438162, 
	     11.2723074,      12.1597347,      13.1093397,      14.1243973, 
	     15.2081976,      16.3640213,      17.5950966,      18.9045696, 
	     20.2954617,      21.7706108,      23.3326225,      24.9838047, 
	     26.7260971,      28.5609989,      30.4894753,       32.511879, 
	     34.6278381,       36.836174,      39.1347733,      41.5204849, 
	     43.9890137,      46.5348015,      49.1509132,       51.828949, 
	     54.5589333,      57.3292618,      60.1266174,       62.935955, 
	     65.7405014,      68.5217972,      71.2597733,      73.9328842, 
	     76.5183411,      78.9923477,      81.3304749,      83.5080414, 
	     85.5006638,      87.2848434,      88.8385849,      90.1422424, 
	     91.1792221,      91.9369507,      92.4076996,      92.5894775, 
	     92.4869461,       92.112114,      91.4850693,      90.6343994, 
	     89.5974579,      88.4202957,      87.1572495,      85.8701706, 
	     84.6270905,      83.5005417,      82.5653152,      81.8957672, 
	     81.5626755,      81.6298218,      82.1502457,      83.1625443, 
	     84.6872253,      86.7234268,      89.2463074,      92.2052536, 
	     95.5233002,      99.0979156,      102.803383,      106.494949, 
	     110.014687,      113.199059,      115.887917,       117.93454, 
	     119.216202,      119.644577,       119.17514,      117.814651, 
	     115.625893,      112.728798,      109.297203,      105.550888, 
	     101.742867,      98.1422424,      95.0134811,      92.5936661, 
	     91.0694427,      90.5558701,      91.0796356,      92.5688553, 
	      94.851326,      97.6625671,      100.663872,       103.46962, 
	     105.681839,      106.929024,      106.905228,       105.40519, 
	     102.351059,      97.8071594,      91.9801865,      85.2037354, 
	     77.9081345,       70.578186,      63.7033157,      57.7258224, 
	     52.9934425,       49.722332,      47.9752312,      47.6580467, 
	     48.5350647,      50.2607956,       52.423893,      54.5968742, 
	     56.3847847,      57.4664154,      57.6231003,       56.752449, 
	      54.866703,      52.0777321,      48.5722427,      44.5817566, 
	     40.3518639,      36.1146812,      32.0673256,      28.3580551, 
	     25.0802326,       22.273489,      19.9305019,      18.0075665, 
	     16.4369278,      15.1391201,      14.0339012,      13.0488214, 
	     12.1250029,      11.2201004,      10.3087597,      9.38113308, 
	     8.44005108,      7.49750996,      6.57099962,      5.68010855, 
	     4.84368563,       4.0777216,      3.39398646,       2.7993803, 
	     2.29589605,      1.88105953,      1.54869676,      1.28988755, 
	     1.09397805,     0.949547112,     0.845253229,     0.770512223, 
	    0.715984941,     0.673874974,     0.638054609,     0.604046702, 
	    0.568897545,     0.530975521,     0.489727706,      0.44542101, 
	    0.398890644,     0.351309568,     0.303989798,     0.258221358, 
	    0.215150192,     0.175695851,     0.140505522,     0.109940849, 
	   0.0840920061,    0.0628127009,    0.0457692482,    0.0324965976, 
	   0.0224549007,     0.015081034,   0.00983107835,   0.00621139724, 
	  0.00379776093,   0.00224342081,   0.00127818808,  0.000701125304, 
	 0.000369562156,  0.000186810314,  9.03692198e-05, 
/* Z=68 -> Er */
	              0,     0.247885138,     0.270068645,     0.294220626, 
	    0.320513338,      0.34913373,     0.380284667,     0.414186209, 
	    0.451077074,     0.491216034,     0.534883559,     0.582383573, 
	    0.634045243,     0.690224946,     0.751308382,     0.817712843, 
	    0.889889479,     0.968326032,      1.05354917,      1.14612806, 
	     1.24667585,      1.35585523,      1.47437906,      1.60301614, 
	     1.74259317,      1.89399993,      2.05819201,       2.2361958, 
	     2.42911196,      2.63812065,      2.86448479,      3.10955477, 
	      3.3747735,      3.66167998,      3.97191358,      4.30721855, 
	     4.66944695,      5.06056309,      5.48264551,      5.93789005, 
	     6.42861128,      6.95724249,      7.52633667,      8.13856411, 
	     8.79671001,      9.50366974,      10.2624388,      11.0761089, 
	      11.947855,      12.8809156,       13.878582,        14.94417, 
	     16.0809975,      17.2923565,      18.5814705,      19.9514599, 
	     21.4052868,       22.945715,      24.5752392,      26.2960224, 
	     28.1098232,      30.0179157,      32.0209961,      34.1191025, 
	     36.3114891,      38.5965538,      40.9717026,      43.4332504, 
	     45.9762917,      48.5946121,       51.280571,      54.0249825, 
	     56.8170624,      59.6443405,      62.4926109,      65.3459167, 
	     68.1865692,      70.9952087,      73.7508926,      76.4313049, 
	     79.0129318,       81.471405,      83.7818604,      85.9194031, 
	     87.8596573,      89.5793915,      91.0572357,       92.274437, 
	     93.2157516,      93.8702927,      94.2324448,      94.3027496, 
	      94.088768,      93.6057968,      92.8774567,      91.9361496, 
	     90.8231125,       89.588295,      88.2897644,      86.9927292, 
	     85.7680969,      84.6904907,      83.8358688,      83.2786102, 
	     83.0881958,      83.3256302,        84.03965,      85.2629471, 
	     87.0085983,      89.2669678,      92.0033264,      95.1564331, 
	     98.6384735,      102.336372,      106.114899,       109.82148, 
	     113.292809,      116.363098,      118.873619,      120.683212, 
	     121.679024,      121.786774,      120.979752,      119.285522, 
	     116.789429,      113.634239,      110.015099,      106.169647, 
	     102.363358,      98.8706131,      95.9526901,      93.8341904, 
	     92.6800308,      92.5752563,      93.5101776,      95.3729477, 
	      97.951416,      100.945183,      103.987663,      106.677094, 
	     108.613876,      109.440826,      108.882103,      106.776329, 
	      103.09967,      97.9754181,      91.6679535,      84.5608902, 
	     77.1208954,       69.850708,      63.2365303,      57.6959839, 
	     53.5330124,       50.905735,      49.8113556,      50.0905151, 
	     51.4502449,      53.5023613,      55.8117828,      57.9479752, 
	     59.5324745,      60.2764397,      60.0040359,      58.6599274, 
	     56.3015747,       53.079174,      49.2074814,      44.9341583, 
	     40.5092125,      36.1590805,      32.0676537,      28.3654728, 
	     25.1266823,      22.3727589,      20.0812244,        18.19734, 
	     16.6468334,      15.3479767,      14.2217922,      13.1996336, 
	     12.2278948,      11.2699728,      10.3059244,      9.33039761, 
	     8.34949017,       7.3771286,       6.4314723,      5.53171396, 
	     4.69550705,       3.9371295,      3.26639199,      2.68821645, 
	     2.20277762,      1.80606163,      1.49069953,      1.24694085, 
	      1.0636462,     0.929210782,     0.832350791,      0.76271534, 
	    0.711310506,     0.670740247,     0.635285735,     0.600853443, 
	     0.56482482,     0.525842369,     0.483561814,     0.438395709, 
	    0.391267449,     0.343389571,     0.296074867,     0.250585109, 
	    0.208018869,     0.169237435,     0.134826481,     0.105088808, 
	   0.0800629705,    0.0595611595,    0.0432197452,    0.0305554178, 
	   0.0210207235,     0.014053788,   0.00911851227,   0.00573327579, 
	  0.00348784449,   0.00204963633,   0.00116148335,  0.000633543765, 
	 0.000331998803,  0.000166808371,  8.01862261e-05, 
/* Z=69 -> Tm */
	              0,      0.26868698,      0.29260397,       0.3186315, 
	    0.346953183,     0.377768129,     0.411292344,     0.447760016, 
	    0.487425029,     0.530562401,     0.577470183,     0.628471017, 
	    0.683914244,     0.744177759,     0.809670389,     0.880834222, 
	    0.958146811,      1.04212415,      1.13332331,       1.2323451, 
	     1.33983791,      1.45650005,      1.58308375,       1.7203989, 
	     1.86931646,      2.03077197,      2.20577049,      2.39539003, 
	     2.60078645,      2.82319713,      3.06394529,      3.32444501, 
	     3.60620546,      3.91083455,      4.24004412,       4.5956521, 
	     4.97958851,      5.39389658,      5.84073591,      6.32238531, 
	     6.84124327,      7.39982796,      8.00077629,      8.64684391, 
	     9.34089661,      10.0859089,      10.8849554,      11.7411947, 
	     12.6578674,      13.6382666,      14.6857309,      15.8036098, 
	     16.9952393,      18.2639122,      19.6128368,      21.0450878, 
	     22.5635681,      24.1709423,      25.8695793,       27.661478, 
	     29.5481949,      31.5307522,      33.6095581,      35.7843056, 
	     38.0538559,      40.4161491,      42.8680878,       45.405407, 
	      48.022583,      50.7127151,      53.4674072,      56.2766953, 
	     59.1289444,      62.0108109,      64.9071808,      67.8011856, 
	     70.6742325,      73.5060577,      76.2749023,      78.9576645, 
	     81.5301895,      83.9675903,      86.2446594,      88.3363647, 
	     90.2184601,      91.8681183,      93.2647018,      94.3905563, 
	     95.2318954,      95.7797318,      96.0307083,      95.9880829, 
	     95.6624603,      95.0725479,      94.2456436,      93.2179871, 
	     92.0347519,      90.7497253,      89.4246063,      88.1278076, 
	     86.9328079,      85.9160385,      85.1541824,      84.7211533, 
	     84.6845932,      85.1021423,       86.017601,      87.4571609, 
	     89.4259109,      91.9049683,      94.8493423,      98.1870117, 
	      101.81929,      105.622871,      109.453491,      113.151588, 
	     116.549553,      119.480621,      121.789032,       123.34079, 
	     124.034523,      123.811562,      122.664307,      120.642021, 
	     117.852982,      114.462547,      110.686211,      106.777763, 
	     103.012703,      99.6676865,      96.9973679,      95.2103348, 
	     94.4464111,      94.7578659,       96.096756,      98.3106842, 
	      101.14846,      104.276199,      107.303459,      109.817574, 
	      111.42337,      111.784264,      110.660347,      107.938828, 
	     103.652748,       97.984848,      91.2551727,       83.892746, 
	     76.3937759,      69.2706146,      62.9973297,      57.9584198, 
	     54.4071732,      52.4392624,      51.9850502,      52.8218842, 
	     54.6044998,       56.909214,      59.2856445,      61.3086548, 
	     62.6236115,      62.9794159,       62.245945,      60.4152222, 
	      57.587944,      53.9489937,      49.7365532,      45.2096519, 
	     40.6183739,      36.1800804,      32.0632782,      28.3798542, 
	     25.1849194,      22.4828434,      20.2375565,      18.3850765, 
	     16.8463192,      15.5387192,       14.385581,      13.3226328, 
	     12.3016949,      11.2917233,      10.2777634,      9.25841713, 
	     8.24247169,      7.24525118,      6.28514433,      5.38063002, 
	     4.54797173,      3.79965496,      3.14354515,      2.58267069, 
	     2.11551547,      1.73667502,      1.43773735,      1.20825922, 
	     1.03673339,     0.911461413,     0.821280301,     0.756110966, 
	    0.707321525,     0.667916536,     0.632575512,     0.597570121, 
	     0.56059581,      0.52054745,      0.47726965,     0.431303084, 
	    0.383645415,     0.335538208,     0.288287699,     0.243123055, 
	    0.201093271,     0.163001314,     0.129372582,     0.100453302, 
	   0.0762331113,    0.0564856417,    0.0408200845,    0.0287372842, 
	   0.0196840856,    0.0131011894,   0.00846108329,   0.00529443799, 
	   0.0032049031,   0.00187368365,   0.00105611281,  0.000572878867, 
	  0.00029848158,  0.000149071348,  7.12141773e-05, 
/* Z=70 -> Yb */
	              0,     0.291215748,     0.316996723,     0.345039815, 
	    0.375540614,     0.408711165,     0.444781303,     0.484000146, 
	    0.526637495,     0.572985709,     0.623360991,      0.67810595, 
	    0.737590909,     0.802216589,     0.872416019,     0.948657155, 
	     1.03144526,      1.12132597,      1.21888769,      1.32476532, 
	     1.43964243,      1.56425583,      1.69939816,      1.84592211, 
	     2.00474405,       2.1768477,      2.36328912,      2.56519961, 
	     2.78379178,       3.0203619,       3.2762959,      3.55307293, 
	     3.85226989,      4.17556572,      4.52474499,      4.90170288, 
	     5.30844641,       5.7470994,      6.21990395,      6.72922182, 
	     7.27753544,      7.86744833,      8.50168133,      9.18307114, 
	     9.91456699,      10.6992178,      11.5401716,      12.4406567, 
	     13.4039688,      14.4334555,      15.5324898,      16.7044487, 
	     17.9526806,      19.2804642,      20.6909828,      22.1872559, 
	     23.7721024,      25.4480801,      27.2174072,      29.0819035, 
	     31.0428944,      33.1011391,      35.2567253,      37.5089722, 
	     39.8563194,      42.2962303,      44.8250656,      47.4379616, 
	     50.1287422,       52.889782,      55.7119293,      58.5843925, 
	     61.4946861,       64.428566,      67.3700104,       70.301239, 
	     73.2027435,      76.0534058,      78.8306427,       81.510643, 
	     84.0686569,      86.4793243,      88.7172012,      90.7572327, 
	      92.575386,      94.1493912,      95.4594803,      96.4892654, 
	     97.2265854,      97.6644974,       97.802124,      97.6455383, 
	       97.20858,      96.5134659,      95.5912704,      94.4821243, 
	     93.2351303,      91.9078751,      90.5655289,      89.2795029, 
	     88.1256027,      87.1816406,      86.5246353,      86.2275085, 
	      86.355484,      86.9622574,      88.0860748,      89.7460403, 
	      91.938736,      94.6355667,      97.7810135,      101.292137, 
	      105.05954,      108.949951,      112.810783,      116.476334, 
	     119.775841,      122.542999,      124.626534,      125.901222, 
	     126.278816,       125.71769,      124.230499,      121.888924, 
	     118.824371,      115.224274,      111.323242,      107.389236, 
	      103.70517,      100.546799,      98.1586838,       96.729866, 
	     96.3719711,      97.1019974,      98.8323441,      101.370018, 
	     104.426247,      107.636703,       110.59153,      112.872742, 
	     114.095741,      113.950691,      112.239052,      108.900612, 
	     104.027222,      97.8606339,      90.7734985,      83.2347717, 
	      75.762764,      68.8708344,      63.0121422,      58.5302277, 
	     55.6218033,      54.3169594,      54.4792366,      55.8259277, 
	     57.9653816,      60.4462013,      62.8111801,      64.6486206, 
	     65.6342697,      65.5590973,      64.3404922,      62.0171585, 
	     58.7303276,      54.6954842,      50.1695328,      45.4182587, 
	     40.6879196,      36.1839409,      32.0577927,      28.4022598, 
	     25.2539482,      22.6013641,      20.3964958,      18.5678291, 
	     17.0330429,      15.7099609,       14.524991,      13.4186268, 
	     12.3481302,      11.2877493,       10.227066,      9.16810036, 
	     8.12178707,      7.10437298,      6.13410282,       5.2284832, 
	     4.40224504,      3.66604018,       3.0258286,      2.48284197, 
	     2.03400445,       1.6726644,      1.38950658,      1.17352092, 
	     1.01293409,     0.896031916,      0.81182456,     0.750533879, 
	    0.703901172,     0.665328562,     0.629880846,     0.594177008, 
	    0.556205809,      0.51509434,     0.470857531,      0.42414856, 
	     0.37602675,     0.327753395,     0.280621678,     0.235824406, 
	    0.194359347,     0.156971171,     0.124126367,    0.0960166901, 
	   0.0725855231,    0.0535705909,    0.0385565162,    0.0270304903, 
	   0.0184353702,    0.0122156125,   0.00785296038,   0.00489058346, 
	  0.00294587901,   0.00171347032,  0.000960698759,  0.000518260116, 
	 0.000268482981,   0.00013329339,  6.32836382e-05, 
/* Z=71 -> Lu */
	              0,     0.315639436,     0.343426526,     0.373637587, 
	    0.406480938,     0.442182213,      0.48098594,     0.523156881, 
	    0.568981826,     0.618771136,     0.672860444,     0.731613159, 
	    0.795421839,     0.864710987,     0.939939141,      1.02160132, 
	     1.11023188,      1.20640755,      1.31074953,      1.42392766, 
	     1.54666328,      1.67973268,      1.82397056,      1.98027432, 
	     2.14960718,      2.33300328,      2.53157091,      2.74649739, 
	     2.97905278,      3.23059511,      3.50257421,      3.79653597, 
	     4.11412716,      4.45709848,      4.82731056,      5.22673464, 
	     5.65745831,      6.12168646,      6.62174368,      7.16007566, 
	     7.73924971,      8.36195183,      9.03098869,      9.74927711, 
	     10.5198441,      11.3458157,      12.2304068,      13.1769104, 
	     14.1886749,      15.2690926,      16.4215698,      17.6494999, 
	     18.9562283,      20.3450184,       21.819006,      23.3811474, 
	     25.0341644,      26.7804832,      28.6221581,      30.5608006, 
	      32.597496,      34.7327003,      36.9661674,      39.2968102, 
	     41.7226257,      44.2405624,      46.8464127,       49.534687, 
	     52.2985229,      55.1295586,      58.0178413,      60.9517479, 
	     63.9179115,      66.9011841,      69.8846283,      72.8495483, 
	     75.7755432,      78.6406479,      81.4215012,      84.0936203, 
	     86.6317139,      89.0100555,      91.2030334,      93.1856613, 
	     94.9342957,      96.4273224,      97.6460114,      98.5753937, 
	     99.2051315,      99.5304794,      99.5531998,      99.2823792, 
	     98.7352066,       97.937561,      96.9243546,      95.7396698, 
	     94.4365234,      93.0762024,      91.7272263,      90.4637909, 
	     89.3636932,      88.5057831,      87.9669113,      87.8185043, 
	     88.1227875,      88.9288788,      90.2688675,      92.1541595, 
	     94.5723038,      97.4845734,      100.824638,      104.498573, 
	      108.38633,      112.345169,      116.214767,      119.824272, 
	     123.000977,      125.580437,       127.41748,      128.397354, 
	     128.446533,      127.541931,      125.717644,      123.068481, 
	     119.749184,      115.968719,      111.979553,      108.061821, 
	     104.503098,      101.575043,      99.5084915,      98.4693298, 
	       98.53759,      99.6923294,      101.804588,      104.640343, 
	     107.874268,      111.114227,       113.93499,      115.918411, 
	     116.696404,      115.991989,      113.653587,      109.678078, 
	     104.218964,      97.5776825,      90.1778488,      82.5243759, 
	     75.1517105,      68.5668106,      63.1938133,      59.3273506, 
	     57.1007118,      56.4733124,      57.2391586,      59.0552483, 
	     61.4857941,      64.0559082,       66.307251,      67.8479919, 
	     68.3908081,      67.7749023,      65.9705048,      63.0674019, 
	     59.2506714,      54.7686691,      49.8982697,      44.9122009, 
	     40.0520935,      35.5096397,      31.4165344,      27.8429279, 
	     24.8028183,      22.2644882,      20.1637325,      18.4177742, 
	     16.9381561,      15.6413288,      14.4562769,      13.3289671, 
	     12.2238798,      11.1231127,      10.0237589,       8.9342432, 
	     7.87029219,      6.85104322,      5.89569187,      5.02089453, 
	     4.23902893,      3.55728388,      2.97748804,        2.496526, 
	     2.10716867,      1.79914939,        1.560323,      1.37778175, 
	     1.23882985,      1.13175476,      1.04636884,      0.97432375, 
	    0.909217775,     0.846535802,     0.783465147,     0.718630672, 
	    0.651791036,     0.583527088,     0.514948726,     0.447436363, 
	    0.382426739,     0.321247786,     0.265001923,     0.214495316, 
	    0.170207649,     0.132296324,      0.10062743,    0.0748258457, 
	   0.0543366186,    0.0384901278,    0.0265645254,    0.0178402215, 
	   0.0116429236,   0.00737349736,   0.00452467892,   0.00268611009, 
	  0.00154015969,  0.000851454271,  0.000453023385,  0.000231533966, 
	 0.000113442118,  5.31726437e-05,  2.37904005e-05, 
/* Z=72 -> Hf */
	              0,     0.342106819,     0.372051656,     0.404593199, 
	    0.439953089,     0.478371471,     0.520108223,     0.565444827, 
	    0.614685714,     0.668160319,      0.72622484,     0.789264262, 
	    0.857694745,     0.931965649,      1.01256216,      1.10000765, 
	     1.19486725,      1.29774952,      1.40931046,      1.53025651, 
	     1.66134775,      1.80340159,      1.95729685,      2.12397671, 
	     2.30445409,      2.49981403,      2.71121955,      2.93991494, 
	     3.18723059,      3.45458746,      3.74350071,      4.05558538, 
	     4.39255857,      4.75624609,      5.14858389,      5.57162237, 
	     6.02752876,      6.51859045,      7.04721546,      7.61593342, 
	     8.22739601,      8.88437176,      9.58975029,      10.3465281, 
	     11.1578083,      12.0267878,       12.956749,      13.9510403, 
	     15.0130606,      16.1462402,      17.3540058,      18.6397591, 
	     20.0068378,      21.4584713,      22.9977417,      24.6275196, 
	     26.3504143,      28.1687031,      30.0842552,      32.0984573, 
	     34.2121201,      36.4253883,      38.7376366,      41.1473618, 
	     43.6520767,       46.248188,      48.9308968,      51.6940575, 
	     54.5300903,      57.4298706,      60.3826332,      63.3758888, 
	     66.3953857,      69.4250565,       72.447052,      75.4417572, 
	     78.3879013,      81.2626877,      84.0420456,      86.7008514, 
	     89.2133408,      91.5535126,      93.6956635,      95.6150284, 
	     97.2884445,      98.6951218,      99.8175278,       100.64225, 
	     101.160957,      101.371277,      101.277771,      100.892693, 
	      100.23671,      99.3394547,      98.2397461,      96.9856339, 
	     95.6339798,      94.2496948,      92.9044571,      91.6749649, 
	     90.6406784,      89.8810577,      89.4723053,      89.4838028, 
	     89.9741745,      90.9873657,      92.5487061,      94.6613617, 
	     97.3033829,      100.425583,      103.950661,      107.773697, 
	     111.764381,      115.770973,      119.626266,       123.15535, 
	     126.184914,      128.553909,      130.124847,      130.795151, 
	     130.507523,      129.258636,      127.104988,      124.164879, 
	     120.616013,      116.687897,      112.649078,      108.789352, 
	      105.39785,      102.738579,      101.024986,      100.396362, 
	     100.898239,      102.469749,      104.939903,      108.034477, 
	     111.394096,      114.602814,      117.225304,      118.849518, 
	     119.130386,      117.829857,      114.848282,      110.242767, 
	     104.229317,      97.1672516,      89.5266647,      81.8416977, 
	     74.6544647,      68.4563522,      63.6335602,      60.4243774, 
	     58.8935623,      58.9278412,      60.2529526,      62.4700203, 
	     65.1060028,      67.6709671,      69.7143021,      70.8721924, 
	     70.9008179,      69.6919556,      67.2706604,      63.7774239, 
	     59.4389725,      54.5329742,      49.3520203,      44.1714478, 
	     39.2242775,      34.6851997,      30.6637878,      27.2061577, 
	     24.3032341,      21.9033966,      19.9271393,      18.2816887, 
	     16.8738422,      15.6199789,      14.4526777,      13.3239813, 
	     12.2056732,      11.0872526,      9.97238541,      8.87457275, 
	     7.81274891,      6.80730057,      5.87685299,      5.03600168, 
	     4.29401922,       3.6544385,      3.11536598,      2.67032003, 
	     2.30938411,      2.02048302,      1.79061985,      1.60695696, 
	     1.45766306,      1.33249843,      1.22313964,      1.12327647, 
	     1.02852345,     0.936204195,     0.845054686,     0.754894674, 
	    0.666300356,     0.580305398,     0.498144716,     0.421050519, 
	    0.350100815,     0.286119401,     0.229621574,     0.180799127, 
	    0.139537036,     0.105454028,    0.0779587403,    0.0563139804, 
	   0.0397019461,    0.0272845943,    0.0182544924,    0.0118733821, 
	  0.00749731483,   0.00458886102,   0.00271818554,   0.00155560672, 
	 0.000858619926,  0.000456223468,  0.000232907667,   0.00011400666, 
	 5.33927923e-05,  2.38704215e-05,  1.01632049e-05, 
/* Z=73 -> Ta */
	              0,     0.370798081,     0.403064162,     0.438111186, 
	     0.47617504,     0.517510951,     0.562395215,     0.611126721, 
	    0.664028645,     0.721450508,     0.783769965,     0.851394951, 
	    0.924766064,      1.00435889,      1.09068656,      1.18430245, 
	     1.28580296,      1.39583087,      1.51507783,      1.64428854, 
	     1.78426349,      1.93586326,      2.10001159,      2.27770019, 
	     2.46999216,      2.67802644,      2.90302181,      3.14628172, 
	     3.40919828,      3.69325709,      4.00004101,      4.33123446, 
	      4.6886282,      5.07412291,      5.48973179,      5.93758488, 
	     6.41993141,      6.93914032,      7.49770403,      8.09823608, 
	     8.74347115,      9.43626118,       10.179575,      10.9764872, 
	     11.8301754,      12.7439032,       13.721015,      14.7649107, 
	      15.879034,      17.0668411,      18.3317757,      19.6772385, 
	     21.1065426,      22.6228733,      24.2292442,      25.9284248, 
	     27.7228928,      29.6147594,      31.6056881,      33.6968155, 
	     35.8886604,       38.181015,      40.5728645,       43.062252, 
	     45.6461754,       48.320488,      51.0797501,      53.9171448, 
	     56.8243446,      59.7914391,      62.8068161,      65.8571243, 
	     68.9271927,      72.0000458,      75.0569153,      78.0772858, 
	     81.0390244,      83.9185638,       86.691124,      89.3310242, 
	     91.8121033,      94.1081619,      96.1935425,      98.0437775, 
	     99.6363373,      100.951424,      101.972832,      102.688889, 
	      103.09343,      103.186638,      102.976021,      102.477127, 
	     101.714256,      100.720871,      99.5397263,      98.2228317, 
	     96.8308563,      95.4321518,      94.1013565,      92.9173965, 
	     91.9610291,      91.3118515,      91.0449371,      91.2270279, 
	     91.9125824,      93.1397476,      94.9265213,      97.2673111, 
	     100.130249,      103.455406,      107.154411,      111.111519, 
	     115.186462,      119.219193,      123.036552,      126.460716, 
	     129.319183,      131.455887,      132.742645,      133.090622, 
	     132.460236,       130.86911,      128.396744,      125.185097, 
	     121.434311,       117.39315,       113.34417,      109.584084, 
	     106.400337,       104.04567,      102.712563,      102.510361, 
	     103.447571,      105.421898,      108.220169,      111.529427, 
	     114.959251,      118.074417,      120.435303,      121.642281, 
	      121.37973,      119.454277,      115.822571,      110.604263, 
	     104.077629,      96.6570206,      88.8535843,      81.2230682, 
	     74.3063583,      68.5699158,       64.353508,      61.8323364, 
	     60.9986458,      61.6657448,      63.4937439,      66.0335464, 
	     68.7827225,      71.2453308,      72.9874725,      73.6813049, 
	     73.1325226,      71.2890701,      68.2319336,      64.1511688, 
	     59.3121376,       54.017662,      48.5713272,      43.2463799, 
	     38.2630005,      33.7754784,      29.8688526,      26.5638332, 
	     23.8278618,      21.5897884,      19.7558022,      18.2244473, 
	     16.8992043,      15.6976576,      14.5569592,      13.4357576, 
	     12.3132029,      11.1858025,      10.0630217,      8.96245384, 
	     7.90519571,      6.91195774,      6.00014448,      5.18203449, 
	      4.4639883,      3.84652948,      3.32508039,      2.89111018, 
	     2.53347039,      2.23972011,      1.99729955,      1.79446208, 
	      1.6209234,      1.46823204,       1.3298912,      1.20128703, 
	     1.07947695,     0.962898791,     0.851045012,     0.744142294, 
	    0.642860591,     0.548067451,     0.460634798,     0.381297737, 
	    0.310563475,     0.248663232,     0.195541084,     0.150871426, 
	     0.11409805,    0.0844872296,    0.0611873493,    0.0432893857, 
	    0.029882418,    0.0201002937,    0.0131568136,   0.00836835802, 
	  0.00516439183,   0.00308747427,   0.00178514805,  0.000996500603, 
	 0.000536069623,  0.000277378887,  0.000137771945,  6.55489785e-05, 
	 2.98074829e-05,  1.29249129e-05,  5.33101593e-06, 
/* Z=74 -> W */
	              0,     0.401910275,      0.43667385,     0.474414945, 
	    0.515384376,     0.559853435,     0.608115673,     0.660488367, 
	    0.717314482,      0.77896452,     0.845838726,     0.918369114, 
	    0.997022092,      1.08230054,      1.17474711,      1.27494597, 
	     1.38352704,      1.50116813,      1.62859845,      1.76660216, 
	      1.9160223,       2.0777638,      2.25279737,       2.4421649, 
	     2.64698124,       2.8684411,      3.10782099,      3.36648488, 
	     3.64588881,      3.94758439,       4.2732234,      4.62456274, 
	     5.00346661,      5.41191244,      5.85199213,      6.32591581, 
	     6.83601522,      7.38474226,      7.97467279,       8.6085043, 
	     9.28905296,      10.0192537,      10.8021526,      11.6409006, 
	      12.538744,       13.499012,      14.5251036,      15.6204681, 
	     16.7885799,      18.0329189,      19.3569355,      20.7640152, 
	     22.2574387,      23.8403339,      25.5156231,      27.2859669, 
	     29.1536884,      31.1207142,      33.1884766,      35.3578415, 
	     37.6290054,         40.0014,      42.4735718,       45.043087, 
	     47.7064171,      50.4588013,      53.2941628,      56.2049599, 
	     59.1821098,      62.2148819,      65.2908096,      68.3956451, 
	     71.5133133,      74.6259079,      77.7137451,      80.7554398, 
	     83.7280121,      86.6071548,      89.3674469,      91.9827042, 
	     94.4264526,      96.6723862,      98.6949997,      100.470268, 
	     101.976418,      103.194794,       104.11068,      104.714348, 
	     105.001923,      104.976318,      104.648148,      104.036423, 
	     103.169121,      102.083611,       100.82666,      99.4541702, 
	     98.0305176,      96.6273727,      95.3220596,      94.1954117, 
	     93.3290787,      92.8023758,      92.6886597,      93.0514908, 
	     93.9405518,      95.3875885,      97.4027176,      99.9710999, 
	     103.050568,      106.570251,      110.430687,      114.505508, 
	     118.644958,       122.68145,      126.436852,      129.731705, 
	     132.395782,      134.279541,      135.265823,      135.281036, 
	     134.304657,      132.376266,      129.598892,      126.138031, 
	     122.215446,      118.097603,      114.078773,      110.459648, 
	     107.522591,       105.50544,       104.57618,       104.81118, 
	     106.179672,      108.536919,      111.627998,      115.103241, 
	     118.544807,      121.503258,      123.540642,      124.276283, 
	     123.430267,      120.859238,      116.579773,      110.775452, 
	      103.78566,      96.0758438,      88.1918182,      80.7026825, 
	     74.1388092,      68.9325104,      65.3690643,      63.5547867, 
	     63.4066238,      64.6651077,      66.9286957,      69.7048416, 
	     72.4703827,      74.7327499,      76.0837479,      76.2389755, 
	     75.0589371,      72.5505219,      68.8508224,      64.1973953, 
	     58.8905678,      53.2535782,      47.5958519,      42.1835861, 
	     37.2198334,      32.8350868,      29.0875473,      25.9713192, 
	     23.4301357,      21.3739338,      19.6957569,      18.2869396, 
	     17.0491467,      15.9024916,      14.7896729,       13.676487, 
	     12.5495396,      11.4120607,      10.2787828,      9.17074299, 
	     8.11062717,       7.1190958,      6.21227026,      5.40040255, 
	     4.68757582,       4.0722146,       3.5481391,      3.10589862, 
	     2.73416018,      2.42097521,      2.15481377,      1.92531025, 
	     1.72371852,      1.54310036,      1.37830675,      1.22580636, 
	     1.08342516,     0.950045168,     0.825301707,     0.709308445, 
	    0.602423251,     0.505062997,      0.41756767,       0.3401106, 
	     0.27264896,     0.214907661,     0.166390285,     0.126409784, 
	   0.0941327214,    0.0686305389,    0.0489326194,    0.0340758637, 
	   0.0231473017,    0.0153168663,   0.00985912606,   0.00616395287, 
	  0.00373727456,   0.00219389284,   0.00124478969,  0.000681421894, 
	   0.0003592172,  0.000181995347,  8.84349793e-05,  4.11247674e-05, 
	 1.82602253e-05,  7.72309886e-06,  3.10361497e-06, 
/* Z=75 -> Re */
	              0,     0.435660452,     0.473111749,     0.513750196, 
	    0.557842493,     0.605676949,     0.657565176,     0.713843822, 
	    0.774876654,     0.841056347,     0.912806869,     0.990585566, 
	     1.07488585,      1.16623914,      1.26521885,      1.37244165, 
	     1.48857236,      1.61432505,      1.75046849,      1.89782798, 
	     2.05729032,      2.22980642,      2.41639662,      2.61815357, 
	     2.83624673,      3.07192731,      3.32653189,      3.60148692, 
	     3.89831376,      4.21863174,      4.56416321,      4.93673849, 
	     5.33829737,      5.77089453,      6.23670197,      6.73801184, 
	     7.27723694,      7.85691309,      8.47969913,      9.14837551, 
	     9.86584187,      10.6351099,      11.4593039,      12.3416471, 
	       13.28545,      14.2941036,      15.3710585,      16.5197983, 
	     17.7438316,      19.0466461,      20.4316921,      21.9023266, 
	     23.4617844,      25.1131191,      26.8591537,      28.7024136, 
	     30.6450539,      32.6887894,      34.8348045,      37.0836678, 
	     39.4352303,      41.8885193,      44.4416389,      47.0916443, 
	     49.8344383,      52.6646385,      55.5754776,      58.5586815, 
	     61.6043777,      64.7009964,      67.8352051,      70.9918365, 
	     74.1539001,      77.3025742,      80.4172745,      83.4757233, 
	     86.4541779,      89.3276062,      92.0699921,       94.654747, 
	     97.0551453,      99.2448883,      101.198738,       102.89325, 
	     104.307549,      105.424271,       106.23037,      106.718201, 
	     106.886368,      106.740692,      106.295013,      105.571907, 
	      104.60318,      103.430153,      102.103561,      100.683197, 
	     99.2370071,      97.8397598,      96.5712433,      95.5138245, 
	     94.7496414,      94.3572006,      94.4076233,      94.9606934, 
	     96.0606995,      97.7324371,      99.9775696,      102.771568, 
	     106.061653,      109.765862,      113.773682,      117.948418, 
	     122.131439,      126.148392,      129.817307,       132.95842, 
	     135.405273,      137.016479,      137.687515,      137.361526, 
	     136.038345,      133.780289,      130.714279,      127.029022, 
	     122.966843,      118.810081,      114.862244,      111.424896, 
	     108.771828,      107.122353,      106.616562,      107.295067, 
	     109.085983,      111.801537,      115.146042,      118.735649, 
	     122.129272,      124.868584,      126.523521,      126.738838, 
	     125.276299,       122.04731,      117.131241,      110.776161, 
	     103.380318,      95.4554214,      87.5748215,      80.3122253, 
	     74.1781235,      69.5617905,      66.6867218,      65.5856476, 
	     66.0987701,      67.8958282,      70.5188675,      73.4396133, 
	     76.1233902,      78.0905457,      78.9671707,       78.519165, 
	     76.6661911,      73.4757385,      69.1401443,       63.941494, 
	     58.2105255,      52.2855072,      46.4761696,      41.0364838, 
	     36.1480064,      31.9141979,      28.3643684,      25.4652061, 
	     23.1370144,      21.2719421,      19.7515945,      18.4620705, 
	     17.3051758,      16.2052593,      15.1118431,      13.9986591, 
	     12.8600683,      11.7059135,      10.5558186,      9.43377686, 
	     8.36359978,      7.36559296,       6.4545188,      5.63877344, 
	     4.92055464,      4.29673433,       3.7601409,      3.30098653, 
	     2.90822816,      2.57071996,      2.27808094,      2.02125978, 
	     1.79281831,      1.58698773,      1.39955497,      1.22764504, 
	     1.06945014,     0.923950672,      0.79065454,     0.669373274, 
	    0.560042262,      0.46258688,     0.376831502,     0.302445233, 
	    0.238920018,     0.185573503,     0.141571254,     0.105962329, 
	   0.0777227655,    0.0558017157,     0.039166037,    0.0268392954, 
	   0.0179326218,    0.0116657643,   0.00737797702,   0.00452944404, 
	  0.00269484683,   0.00155120948,  0.000862353307,  0.000462139898, 
	 0.000238282417,  0.000117966047,  5.59554646e-05,  2.53733433e-05, 
	 1.09734747e-05,  4.51518599e-06,  1.76300284e-06, 
/* Z=76 -> Os */
	              0,     0.472284257,     0.512628317,     0.556383252, 
	    0.603832424,     0.655281961,      0.71106267,     0.771531701, 
	    0.837074578,     0.908107281,     0.985078752,      1.06847274, 
	     1.15881097,      1.25665522,      1.36261046,      1.47732794, 
	     1.60150802,      1.73590374,      1.88132393,       2.0386374, 
	       2.208776,      2.39273882,      2.59159684,       2.8064959, 
	     3.03866148,      3.28940344,      3.56011987,       3.8523016, 
	     4.16753674,        4.507514,      4.87402821,      5.26898432, 
	      5.6943984,      6.15240479,      6.64525509,      7.17532349, 
	     7.74510622,       8.3572216,      9.01441193,      9.71953964, 
	     10.4755831,       11.285634,      12.1528893,      13.0806398, 
	     14.0722599,      15.1311951,      16.2609367,      17.4650059, 
	     18.7469273,      20.1101952,      21.5582371,      23.0943813, 
	     24.7217941,       26.443449,       28.262043,      30.1799507, 
	     32.1991386,      34.3210907,      36.5467148,      38.8762627, 
	     41.3092079,      43.8441582,      46.4787331,      49.2094498, 
	     52.0316162,      54.9391937,      57.9247093,       60.979126, 
	     64.0917587,      67.2501678,      70.4401398,      73.6455917, 
	     76.8486099,      80.0294571,      83.1666489,      86.2370758, 
	     89.2162323,      92.0784073,      94.7970886,      97.3453217, 
	     99.6962509,      101.823677,      103.702766,      105.310768, 
	     106.627876,      107.638161,      108.330437,      108.699287, 
	     108.745987,      108.479378,       107.91671,      107.084221, 
	     106.017609,      104.762192,      103.372681,      101.912651, 
	     100.453484,      99.0728455,      97.8526535,      96.8764725, 
	     96.2264557,      95.9797974,      96.2048264,       96.956955, 
	     98.2744827,      100.174667,      102.650238,      105.666557, 
	     109.159981,      113.037361,      117.177299,      121.433182, 
	     125.638107,      129.611923,      133.169983,      136.133636, 
	     138.341751,      139.662704,      140.006165,      139.333527, 
	     137.665924,      135.089233,      131.754257,      127.872337, 
	     123.705048,      119.548523,      115.712738,      112.496902, 
	     110.162582,      108.907051,      108.839302,      109.961647, 
	     112.159714,      115.202682,      118.755669,      122.403816, 
	     125.687439,      128.145096,      129.361115,      129.011993, 
	      126.90686,      123.016037,       117.48381,      110.622177, 
	      102.88504,      94.8245392,      87.0336227,      80.0814438, 
	     74.4491196,      70.4743805,       68.312439,      67.9188309, 
	     69.0568542,       71.328804,      74.2266769,      77.1950836, 
	      79.697525,      81.2767334,      81.6014023,       80.493782, 
	     77.9361877,       74.057457,      69.1034012,      63.3968887, 
	     57.2940025,       51.142334,      45.2462425,      39.8423691, 
	     35.0866203,      31.0524673,      27.7386913,      25.0841007, 
	     22.9861641,      21.3205986,      19.9593716,      18.7852993, 
	     17.7021332,      16.6399307,      15.5560427,      14.4326782, 
	     13.2720747,      12.0905342,      10.9122734,      9.76396561, 
	     8.67042351,       7.6516819,      6.72146606,      5.88684177, 
	     5.14877558,      4.50327826,      3.94281864,       3.4577632, 
	     3.03765655,       2.6722374,      2.35214877,       2.0693543, 
	      1.8173033,      1.59090614,      1.38638151,      1.20103586, 
	     1.03301597,     0.881070793,     0.744339824,     0.622179031, 
	    0.514026701,     0.419307739,     0.337372065,     0.267462283, 
	    0.208704874,     0.160120383,     0.120647117,    0.0891741961, 
	   0.0645788237,      0.04576426,    0.0316944607,    0.0214226581, 
	   0.0141118784,   0.00904649962,   0.00563503336,    0.0034051782, 
	  0.00199290877,   0.00112767436,  0.000615799625,  0.000323913235, 
	 0.000163789955,  7.94530788e-05,  3.68936489e-05,  1.63612931e-05, 
	 6.91303558e-06,  2.77597883e-06,  1.05661775e-06, 
/* Z=77 -> Ir */
	              0,     0.512041748,     0.555499911,     0.602607548, 
	    0.653665662,     0.708999336,     0.768959403,     0.833924472, 
	    0.904303193,     0.980536044,      1.06309795,      1.15250075, 
	     1.24929535,      1.35407543,      1.46747935,      1.59019375, 
	      1.7229569,      1.86656165,      2.02185917,      2.18976283, 
	     2.37125206,      2.56737542,      2.77925611,      3.00809479, 
	     3.25517559,      3.52186823,      3.80963421,      4.12002993, 
	     4.45471239,      4.81544161,      5.20408583,      5.62262535, 
	     6.07315493,      6.55788803,      7.07915831,      7.63942099, 
	     8.24125576,      8.88736439,      9.58057117,      10.3238182, 
	     11.1201639,      11.9727736,      12.8849154,      13.8599453, 
	     14.9012947,      16.0124569,      17.1969624,      18.4583569, 
	     19.8001728,      21.2259007,      22.7389374,      24.3425617, 
	     26.0398674,      27.8337154,       29.726675,      31.7209435, 
	     33.8182755,      36.0199051,      38.3264389,      40.7377815, 
	     43.2530098,      45.8702736,       48.586689,       51.398201, 
	     54.2994957,      57.2838516,      60.3430595,      63.4672928, 
	     66.6450195,      69.8629532,      73.1059494,      76.3570099, 
	     79.5973053,      82.8061676,      85.9612579,      89.0386581, 
	     92.0131149,      94.8583374,      97.5473404,      100.052933, 
	     102.348183,      104.407127,       106.20546,      107.721268, 
	     108.935989,      109.835258,      110.409927,      110.656998, 
	     110.580544,      110.192596,      109.513939,      108.574577, 
	     107.414192,      106.082069,      104.636871,      103.145866, 
	     101.683708,      100.330688,      99.1705627,      98.2876587, 
	     97.7636871,       97.673996,      98.0835648,      99.0428162, 
	     100.583481,      102.714752,      105.419899,       108.65387, 
	     112.341927,      116.379753,      120.635277,      124.952454, 
	     129.156876,      133.063583,      136.486557,       139.24968, 
	      141.19873,      142.213486,       142.21936,      141.197205, 
	     139.190598,      136.309235,      132.727829,      128.679428, 
	      124.44326,      120.326973,      116.644058,      113.687927, 
	     111.704376,      110.865143,      111.245071,      112.806122, 
	     115.390419,      118.724701,      122.436974,      126.085129, 
	     129.195923,      131.311005,      132.035446,      131.083817, 
	     128.317795,      123.770103,      117.651077,       110.33519, 
	     102.327759,      94.2146378,      86.5998764,      80.0385513, 
	     74.9731369,      71.6818619,      70.2460861,      70.5418243, 
	     72.2566376,      74.9300079,      78.0113449,      80.9274368, 
	     83.1498871,      84.2531815,      83.9562836,      82.1430817, 
	     78.8610001,      74.2999649,      68.7565689,      62.5899315, 
	     56.1758118,      49.8650627,      43.9508133,      38.6472816, 
	     34.0807686,      30.2921028,      27.2482586,       24.860157, 
	     23.0034275,      21.5390759,      20.3315468,      19.2625256, 
	     18.2396584,      17.2002468,      16.1105442,      14.9618168, 
	     13.7643957,      12.5409946,      11.3202877,      10.1315212, 
	     9.00051975,       7.9472332,      6.98467922,         6.11903, 
	     5.35049248,      4.67463827,      4.08388758,      3.56891561, 
	     3.11983395,      2.72708201,      2.38200855,       2.0771842, 
	     1.80649626,      1.56508744,      1.34920168,      1.15598524, 
	    0.983276308,     0.829411745,     0.693059206,     0.573082864, 
	    0.468441665,     0.378117651,     0.301070571,     0.236214757, 
	    0.182414308,     0.138492405,     0.103251107,    0.0754979253, 
	   0.0540752448,    0.0378893912,    0.0259360895,    0.0173200611, 
	   0.0112672942,   0.00712947734,   0.00438106898,   0.00261020847, 
	  0.00150522531,  0.000838659762,   0.00045063076,  0.000233057566, 
	  0.00011577941,  5.51319026e-05,   2.5108011e-05,  1.09106177e-05, 
	 4.51289225e-06,  1.77223546e-06,  6.59003774e-07, 
/* Z=78 -> Pt */
	              0,     0.555177689,     0.601985514,     0.652697086, 
	    0.707631946,     0.767135024,     0.831578434,     0.901363373, 
	    0.976922512,      1.05872214,      1.14726448,      1.24309051, 
	     1.34678233,      1.45896626,      1.58031559,      1.71155417, 
	     1.85345936,      2.00686574,      2.17266822,      2.35182643, 
	     2.54536867,      2.75439525,      2.98008323,      3.22369051, 
	     3.48655963,        3.770123,      4.07590723,       4.4055357, 
	     4.76073551,      5.14333916,      5.55529022,      5.99864531, 
	     6.47557831,      6.98838282,      7.53947401,      8.13138962, 
	     8.76679134,      9.44846439,      10.1793118,      10.9623556, 
	     11.8007317,      12.6976757,      13.6565237,      14.6806898, 
	     15.7736597,      16.9389668,      18.1801739,      19.5008392, 
	     20.9044991,      22.3946209,      23.9745655,      25.6475449, 
	     27.4165592,      29.2843533,       31.253336,      33.3255081, 
	     35.5024033,      37.7849655,      40.1734886,      42.6674957, 
	     45.2656403,      47.9655914,      50.7639236,      53.6559944, 
	     56.6358299,      59.6960068,      62.8275414,      66.0198059, 
	     69.2604218,      72.5351868,      75.8280563,      79.1211166, 
	     82.3946075,      85.6269684,      88.7950058,      91.8740311, 
	     94.8381195,      97.6604309,      100.313622,      102.770287, 
	     105.003609,      106.987938,      108.699631,      110.117783, 
	     111.225182,      112.009262,      112.463043,      112.586067, 
	     112.385429,       111.87648,      111.083611,      110.040703, 
	     108.791412,      107.389053,      105.896133,      104.383453, 
	     102.928719,      101.614616,      100.526299,      99.7484818, 
	      99.361908,      99.4395294,      100.042374,      101.215302, 
	     102.982918,      105.345802,       108.27742,         111.722, 
	     115.593628,      119.776962,      124.129646,      128.486816, 
	     132.667526,      136.483185,      139.747757,      142.289276, 
	     143.962051,      144.658875,      144.322327,      142.953918, 
	     140.620316,      137.455307,      133.656784,      129.478287, 
	      125.21463,      121.182137,      117.694489,      115.035591, 
	     113.431732,       113.02565,        113.8554,      115.840881, 
	     118.780533,      122.360008,      126.173317,      129.755768, 
	     132.626511,      134.336823,      134.519394,      132.932922, 
	     129.496155,      124.306808,      117.641121,      109.933022, 
	     101.733551,      93.6543884,      86.3012924,      80.2056198, 
	     75.7623901,      73.1835556,      72.4728394,      73.4264374, 
	     75.6595154,      78.6551285,      81.8287125,      84.5988464, 
	     86.4544907,      87.0096664,      86.0389709,      83.4905014, 
	     79.4766159,      74.2457199,      68.1409912,      61.5527878, 
	     54.8716164,      48.4473801,      42.5590096,      37.3967781, 
	     33.0572739,       29.549736,      26.8110809,      24.7263012, 
	     23.1508579,      21.9319706,      20.9264545,      20.0136127, 
	     19.1026878,      18.1351147,      17.0825195,      15.9417381, 
	     14.7281723,       13.468791,      12.1956797,      10.9408369, 
	      9.7324543,      8.59275532,       7.5371418,      6.57438898, 
	     5.70751047,      4.93498421,      4.25206757,      3.65201712, 
	     3.12709999,      2.66935635,      2.27111197,      1.92527878, 
	      1.6254921,      1.36613345,      1.14228725,      0.94966203, 
	    0.784502804,     0.643508136,     0.523757637,     0.422652811, 
	     0.33787024,     0.267323971,     0.209135845,     0.161612287, 
	    0.123225428,    0.0925985277,    0.0684947148,    0.0498088822, 
	   0.0355612263,    0.0248924717,    0.0170592293,    0.0114289857, 
	  0.00747389579,    0.0047630854,   0.00295337732,   0.00177868886, 
	   0.0010386405,  0.000586976879,  0.000320435734,  0.000168642684, 
	 8.53901001e-05,  4.15078503e-05,  1.93270116e-05,  8.60000455e-06, 
	 3.64816219e-06,  1.47159074e-06,  5.62967443e-07, 
/* Z=79 -> Au */
	              0,     0.602082193,     0.652500808,     0.707094669, 
	    0.766202986,     0.830191493,     0.899454713,     0.974417686, 
	     1.05553854,      1.14331043,      1.23826444,      1.34097195, 
	     1.45204759,      1.57215226,      1.70199573,      1.84234035, 
	     1.99400461,        2.157866,      2.33486581,      2.52601147, 
	     2.73238158,      2.95513034,      3.19549036,      3.45477819, 
	     3.73439789,      4.03584576,      4.36071491,      4.71069813, 
	     5.08759403,      5.49330902,       5.9298625,      6.39938974, 
	     6.90414429,      7.44650126,      8.02895641,      8.65413189, 
	     9.32477093,      10.0437384,      10.8140192,      11.6387119, 
	     12.5210238,      13.4642639,      14.4718285,      15.5471859, 
	     16.6938705,      17.9154472,      19.2155037,      20.5976067, 
	      22.065279,      23.6219578,      25.2709503,      27.0153904, 
	     28.8581696,      30.8018913,      32.8487854,      35.0006523, 
	     37.2587509,      39.6237373,      42.0955505,      44.6733093, 
	     47.3552208,      50.1384354,      53.0189629,       55.991539, 
	     59.0494995,      62.1846924,       65.387352,      68.6459961, 
	     71.9473724,      75.2763977,      78.6160812,      81.9475861, 
	     85.2502289,      88.5015869,      91.6776352,      94.7529373, 
	     97.7009506,      100.494362,      103.105499,      105.506927, 
	     107.671982,      109.575478,       111.19455,      112.509483, 
	      113.50457,      114.169197,      114.498711,      114.495438, 
	     114.169548,      113.539886,      112.634598,      111.491501, 
	     110.158272,      108.692207,      107.159622,       105.63472, 
	     104.198059,      102.934357,      101.929878,      101.269241, 
	     101.031738,      101.287369,      102.092628,      103.486229, 
	     105.485062,      108.080597,      111.236076,       114.88472, 
	     118.929359,      123.243683,       127.67543,      132.051514, 
	     136.185303,      139.885681,      142.967957,      145.265732, 
	     146.643433,      147.008408,      146.321793,      144.606903, 
	     141.954163,      138.521561,      134.529816,      130.251724, 
	     125.995781,      122.084518,      118.828743,      116.499512, 
	      115.30024,      115.341621,      116.622482,      119.019257, 
	     122.286652,       126.07058,      129.933945,      133.393829, 
	     135.967163,      137.221008,      136.821716,      134.577469, 
	     130.468277,      124.658791,      117.490761,      109.454536, 
	     101.141312,      93.1808319,      86.1715546,      80.6117172, 
	     76.8403015,      74.9962845,      75.0018997,      76.5730438, 
	     79.2555618,      82.4824448,      85.6439514,       88.160675, 
	     89.5494003,      89.4731827,      87.7698288,      84.4566879, 
	     79.7130737,      73.8450165,      67.2386932,       60.309967, 
	     53.4567032,      47.0196152,      41.2550163,      36.3213272, 
	     32.2786369,      29.0994415,      26.6874275,      24.9006214, 
	     23.5753708,      22.5480499,      21.6722889,      20.8305187, 
	     19.9396038,      18.9512005,      17.8479996,      16.6373444, 
	     15.3436489,      14.0009184,      12.6462107,      11.3146381, 
	     10.0360317,      8.83320141,      7.72150183,      6.70934534, 
	     5.79929972,      4.98943233,      4.27467108,      3.64801121, 
	      3.1014874,      2.62689066,      2.21624494,       1.8620882, 
	     1.55760753,      1.29667664,       1.0738349,     0.884237289, 
	    0.723593891,      0.58810997,     0.474430948,     0.379593521, 
	    0.300982535,     0.236292213,     0.183491185,     0.140790984, 
	    0.106617957,    0.0795888156,    0.0584899113,    0.0422602072, 
	   0.0299774166,    0.0208466426,    0.0141906943,   0.00944109168, 
	  0.00612909719,   0.00387623371,   0.00238409918,   0.00142356567, 
	 0.000823724258,  0.000461024669,  0.000249092293,  0.000129663109, 
	 6.48909845e-05,  3.11543954e-05,  1.43165189e-05,  6.28222915e-06, 
	 2.62590174e-06,  1.04283492e-06,  3.92427381e-07, 
/* Z=80 -> Hg */
	              0,     0.653117001,     0.707429707,     0.766207516, 
	    0.829810381,     0.898626328,     0.973073184,      1.05360126, 
	     1.14069533,      1.23487687,      1.33670723,      1.44678974, 
	     1.56577325,      1.69435441,      1.83328116,      1.98335671, 
	     2.14544177,      2.32045889,      2.50939703,      2.71331358, 
	     2.93333983,      3.17068505,      3.42663956,      3.70258069, 
	     3.99997544,      4.32038641,      4.66547394,      5.03700256, 
	     5.43684387,      5.86698008,      6.32950783,      6.82664156, 
	     7.36071491,      7.93418503,      8.54963112,       9.2097559, 
	     9.91738415,      10.6754627,      11.4870539,      12.3553324, 
	     13.2835741,      14.2751541,      15.3335257,      16.4622116, 
	     17.6647816,      18.9448318,      20.3059597,        21.75173, 
	     23.2856445,      24.9110947,      26.6313248,      28.4493675, 
	     30.3679924,      32.3896446,      34.5163574,      36.7496872, 
	     39.0906181,      41.5394821,      44.0958405,      46.7583809, 
	     49.5248184,      52.3917809,      55.3546638,      58.4075546, 
	     61.5430832,       64.752327,      68.0246964,      71.3478699, 
	      74.707695,      78.0881577,      81.4713516,      84.8375092, 
	     88.1650314,      91.4306335,      94.6095047,      97.6755066, 
	     100.601517,      103.359833,      105.922569,      108.262299, 
	     110.352654,      112.169052,       113.68956,      114.895767, 
	     115.773689,      116.314789,      116.516953,      116.385399, 
	     115.933594,       115.18396,       114.16851,      112.929077, 
	     111.517387,       109.99472,      108.431038,      106.903816, 
	     105.496223,      104.294724,      103.386269,      102.854927, 
	     102.778023,      103.222038,      104.238342,      105.858887, 
	     108.092323,      110.920509,      114.296059,      118.140968, 
	     122.346741,      126.776299,      131.267853,      135.640808, 
	     139.703735,      143.264221,      146.140244,      148.172562, 
	     149.237244,      149.257736,      148.215042,      146.155273, 
	     143.193222,      139.511093,      135.351761,      131.005936, 
	     126.793747,      123.041206,      120.053131,      118.084343, 
	     117.312027,      117.811951,      119.541489,      122.332619, 
	       125.8965,      129.841171,      133.701782,      136.981812, 
	     139.201721,        139.9505,      138.934189,      136.015625, 
	     131.239594,      124.839058,      117.220131,      108.925636, 
	      100.58049,       92.824173,      86.2383575,       81.278656, 
	     78.2196884,      77.1213989,      77.8216248,      79.9567795, 
	     83.0078049,      86.3652878,      89.4043121,       91.558609, 
	     92.3837738,      91.6014786,       89.120018,       85.030098, 
	     79.5788269,      73.1273346,      66.0993652,      58.9284744, 
	     52.0105629,      45.6671562,       40.122261,      35.4938774, 
	     31.7989502,      28.9690418,      26.8731709,      25.3437366, 
	      24.201807,      23.2787304,      22.4321041,      21.5552082, 
	     20.5801258,      19.4755611,      18.2408733,      16.8979912, 
	     15.4827757,       14.037096,      12.6023617,      11.2149696, 
	      9.9036293,      8.68836308,      7.58077049,        6.585145, 
	     5.70001888,      4.91981459,      4.23637438,      3.64023876, 
	     3.12161422,      2.67103839,      2.27978587,      1.94006729, 
	     1.64508784,      1.38900959,       1.1668644,     0.974441111, 
	     0.80816859,     0.665002763,     0.542325199,     0.437853307, 
	     0.34956491,      0.27563706,     0.214398846,     0.164298162, 
	    0.123881347,    0.0917840898,    0.0667314455,    0.0475440696, 
	   0.0331476778,    0.0225826465,    0.0150113637,   0.00972143374, 
	  0.00612388831,   0.00374637847,   0.00222208258,   0.00127562671, 
	 0.000707495841,  0.000378400407,  0.000194788663,   9.6312011e-05, 
	 4.56438647e-05,  2.06875557e-05,  8.94653203e-06,  3.68270366e-06, 
	 1.43927195e-06,  5.32633464e-07,  1.86129313e-07, 
/* Z=81 -> Tl */
	              0,     0.708591223,     0.767097771,     0.830378532, 
	    0.898815334,     0.972819388,      1.05283344,      1.13933384, 
	     1.23283362,      1.33388412,      1.44307792,      1.56105232, 
	     1.68849123,      1.82612908,      1.97475386,      2.13521051, 
	     2.30840421,      2.49530458,      2.69694901,      2.91444731, 
	     3.14898539,      3.40182853,      3.67432761,       3.9679215, 
	     4.28414249,      4.62462044,      4.99108505,      5.38537407, 
	     5.80943346,      6.26532078,      6.75521278,      7.28140259, 
	     7.84630537,      8.45246029,      9.10252762,      9.79929256, 
	       10.54566,       11.344656,       12.199419,      13.1131945, 
	     14.0893269,      15.1312513,      16.2424755,      17.4265614, 
	     18.6871166,      20.0277538,      21.4520721,      22.9636269, 
	     24.5658836,      26.2621746,      28.0556622,      29.9492683, 
	     31.9456177,      34.0469704,      36.2551498,      38.5714493, 
	     40.9965553,      43.5304298,      46.1722412,      48.9202194, 
	     51.7715645,      54.7223358,      57.7673073,      60.8998756, 
	     64.1119308,      67.3937683,       70.733963,      74.1193085, 
	     77.5347519,       80.963356,      84.3862686,      87.7828064, 
	     91.1304855,      94.4051895,       97.581337,      100.632156, 
	     103.530037,      106.246918,      108.754822,      111.026428, 
	     113.035812,      114.759117,      116.175522,      117.268059, 
	      118.02462,      118.438942,      118.511589,      118.250847, 
	     117.673607,      116.806023,      115.683998,      114.353439, 
	     112.870117,      111.299179,      109.714157,      108.195557, 
	     106.828842,      105.701859,       104.90184,      104.511757, 
	     104.606407,      105.248169,      106.482658,      108.334518, 
	     110.803589,      113.861702,      117.450562,      121.480766, 
	     125.832642,      130.358643,      134.888092,      139.233704, 
	     143.200455,      146.595947,      149.242447,      150.989426, 
	     151.726303,      151.394073,      149.994873,      147.598419, 
	     144.344147,      140.438095,      136.144135,      131.769012, 
	     127.641907,      124.089142,      121.405899,      119.827171, 
	     119.500488,      120.463837,      122.631439,      125.790283, 
	     129.609161,       133.66098,      137.457504,      140.494263, 
	     142.301727,      142.497421,      140.833176,      137.231247, 
	     131.803833,      124.852402,      116.844597,      108.370247, 
	     100.080063,      92.6137619,      86.5262604,      82.2212906, 
	     79.9018097,      79.5444031,      80.9019012,      83.5344086, 
	     86.8647385,      90.2504196,      93.0621872,      94.7577896, 
	     94.9412155,      93.4000168,      90.1170502,      85.2572021, 
	     79.1329651,      72.1555405,      64.7792969,      57.4475975, 
	     50.5464325,      44.3709526,       39.107048,      34.8283997, 
	     31.5069218,      29.0333805,      27.2439957,      25.9486713, 
	     24.9570637,      24.0995502,      23.2414398,      22.2899361, 
	     21.1945705,      19.9424343,      18.5500946,      17.0540257, 
	      15.501112,      13.9405022,      12.4173698,      10.9689121, 
	     9.62237549,      8.39478588,      7.29388046,       6.3197546, 
	      5.4668107,      4.72567272,      4.08487177,      3.53219151, 
	      3.0556426,      2.64410472,      2.28768229,      1.97785199, 
	     1.70745575,      1.47060394,      1.26252389,      1.07938707, 
	    0.918132663,     0.776301861,     0.651890218,     0.543221176, 
	    0.448844343,      0.36745742,     0.297850847,     0.238873079, 
	    0.189412564,     0.148392946,     0.114776947,    0.0875757039, 
	   0.0658596382,     0.048768919,    0.0355216749,    0.0254192613, 
	   0.0178483054,    0.0122798849,   0.00826626644,   0.00543569494, 
	  0.00348582701,    0.0021762203,   0.00132024661,  0.000776861503, 
	 0.000442508841,  0.000243509989,  0.000129189561,  6.59362922e-05, 
	  3.2303371e-05,  1.51566601e-05,   6.7945557e-06, 
/* Z=82 -> Pb */
	              0,     0.768951714,     0.831978738,     0.900109828, 
	    0.973749697,      1.05333412,       1.1393317,      1.23224688, 
	     1.33262157,      1.44103873,      1.55812395,      1.68454945, 
	     1.82103634,      1.96835828,      2.12734485,      2.29888439, 
	     2.48392868,      2.68349576,      2.89867401,      3.13062692, 
	     3.38059616,       3.6499064,      3.93996978,      4.25228882, 
	     4.58846235,      4.95018911,      5.33927155,      5.75761938, 
	     6.20725441,      6.69031382,      7.20905113,      7.76584101, 
	     8.36318111,      9.00368977,      9.69011116,      10.4253111, 
	     11.2122717,      12.0540981,      12.9540005,      13.9152937, 
	     14.9413862,      16.0357666,       17.201992,      18.4436626, 
	     19.7644043,      21.1678448,      22.6575775,      24.2371254, 
	     25.9099121,      27.6791992,      29.5480442,      31.5192413, 
	      33.595253,      35.7781372,      38.0694809,      40.4702873, 
	     42.9809074,      45.6009407,      48.3291054,      51.1631584, 
	     54.0997543,      57.1343422,       60.261055,      63.4725723, 
	      66.760025,      70.1128769,      73.5188675,      76.9638824, 
	     80.4319534,      83.9052048,      87.3638687,      90.7863388, 
	     94.1492615,      97.4277115,      100.595413,      103.624992, 
	     106.488411,      109.157349,      111.603844,      113.800819, 
	      115.72287,      117.347061,      118.653793,      119.627724, 
	     120.258797,      120.543198,      120.484322,      120.093681, 
	     119.391708,      118.408386,      117.183617,      115.767403, 
	      114.21949,      112.608833,       111.01236,      109.513359, 
	      108.19931,      107.159058,      106.479607,      106.242371, 
	      106.51902,      107.367233,      108.826286,      110.912933, 
	     113.617706,      116.902023,      120.696404,      124.900032, 
	     129.382126,      133.985153,      138.530136,      142.824127, 
	     146.669601,      149.875671,       152.27034,      153.713501, 
	     154.109467,      153.418167,      151.664017,      148.941025, 
	     145.413315,      141.310196,      136.915253,       132.54921, 
	     128.547546,       125.23365,      122.889587,      121.726891, 
	      121.86026,      123.287453,      125.878082,      129.374207, 
	     133.403793,      137.507843,      141.179398,      143.911865, 
	     145.252213,      144.852966,      142.517395,      138.231232, 
	     132.176117,      124.721184,      116.392097,      107.819366, 
	     99.6710739,      92.5774689,      87.0567627,      83.4520187, 
	     81.8879776,      82.2546234,      84.2202682,      87.2731171, 
	     90.7855911,      94.0923996,      96.5712509,      97.7149734, 
	     97.1854095,      94.8427582,       90.747963,      85.1399994, 
	     78.3933334,      70.9630356,      63.3259583,      55.9260864, 
	     49.1303291,      43.1992035,      38.2739754,      34.3797302, 
	     31.4416809,      29.3106785,      27.7932796,      26.6817837, 
	     25.7803898,      24.9248009,      23.9939938,      22.9142017, 
	     21.6562424,      20.2280445,      18.6644211,      17.0161514, 
	     15.3398867,      13.6900082,       12.112915,      10.6438036, 
	     9.30561924,      8.10967922,      7.05740261,      6.14260006, 
	     5.35388231,      4.67688799,        4.096138,      3.59643579, 
	     3.16382527,      2.78614497,      2.45326853,      2.15710473, 
	     1.89143717,      1.65166974,      1.43452656,      1.23774445, 
	     1.05978346,     0.899573565,     0.756305814,     0.629271746, 
	    0.517750442,     0.420939028,     0.337920159,     0.267658383, 
	    0.209017664,     0.160792246,     0.121743686,    0.0906392708, 
	   0.0662873685,    0.0475674234,    0.0334530808,    0.0230280012, 
	   0.0154945804,    0.0101760952,   0.00651329523,   0.00405645231, 
	  0.00245411089,   0.00143975718,  0.000817609252,  0.000448586885, 
	 0.000237324144,  0.000120822639,  5.90671807e-05,  2.76680057e-05, 
	 1.23892378e-05,  5.29059207e-06,  2.14915758e-06, 
/* Z=83 -> Bi */
	              0,     0.834659696,     0.902560532,     0.975917161, 
	     1.05515862,      1.14074636,      1.23317659,      1.33298254, 
	     1.44073713,      1.55705595,      1.68259919,      1.81807578, 
	      1.9642458,      2.12192369,      2.29198194,      2.47535491, 
	     2.67304182,      2.88611126,      3.11570477,      3.36304092, 
	      3.6294198,      3.91622591,      4.22493505,      4.55711555, 
	     4.91443443,      5.29866028,      5.71166897,       6.1554451, 
	     6.63208771,       7.1438117,       7.6929512,      8.28196144, 
	     8.91342163,      9.59003258,      10.3146172,      11.0901213, 
	     11.9196072,      12.8062496,      13.7533331,      14.7642345, 
	     15.8424244,      16.9914398,      18.2148762,      19.5163631, 
	     20.8995476,      22.3680534,      23.9254532,      25.5752373, 
	     27.3207626,      29.1652031,      31.1115036,      33.1623077, 
	     35.3198891,      37.5860939,      39.9622307,      42.4490013, 
	     45.0464058,      47.7536201,      50.5689163,      53.4895248, 
	     56.5115356,      59.6297684,      62.8376656,      66.1271744, 
	     69.4886246,      72.9106598,      76.3801193,      79.8820038, 
	     83.3994064,      86.9135056,      90.4036102,      93.8472214, 
	     97.2201462,      100.496689,      103.649902,      106.651894, 
	     109.474251,      112.088547,      114.466873,      116.582527, 
	      118.41082,       119.92981,      121.121323,      121.971817, 
	     122.473442,      122.624992,      122.432854,      121.911919, 
	     121.086304,      119.989891,      118.666649,      117.170647, 
	     115.565628,      113.924202,      112.326508,      110.858383, 
	      109.60894,      108.667641,      108.120796,      108.047699, 
	     108.516342,      109.579086,      111.268311,      113.592308, 
	     116.531853,      120.037582,      124.028641,      128.392807, 
	     132.988434,      137.648453,      142.186356,      146.404526, 
	     150.104233,      153.097443,      155.219559,      156.342438, 
	     156.386795,      155.332855,      153.228226,      150.191666, 
	      146.41188,      142.140671,      137.679688,      133.361435, 
	     129.524826,      126.486954,      124.513237,      123.788338, 
	     124.391197,      126.277229,      129.270584,      133.068863, 
	     137.261398,      141.360947,      144.846863,      147.216599, 
	     148.039902,      147.010513,      143.988174,      139.025497, 
	      132.37442,      124.470062,      115.891479,      107.302856, 
	     99.3807602,      92.7362289,      87.8413773,      84.9707108, 
	     84.1655045,      85.2273712,      87.7421799,      91.1319656, 
	     94.7274551,      97.8510895,      99.8986435,      100.408279, 
	     99.1078796,      95.9355087,      91.0319824,      84.7087936, 
	     77.3977356,      69.5907593,      61.7787781,      54.3975372, 
	     47.7864647,      42.1644211,      37.6229706,      34.1359482, 
	     31.5816536,      29.7730942,      28.4910851,      27.5154591, 
	     26.6506233,      25.7431316,      24.6904469,      23.4415455, 
	     21.9909782,      20.3686352,      18.6275444,      16.8318062, 
	     15.0461512,      13.3280897,      11.7229414,      10.2616034, 
	     8.96056747,      7.82355976,      6.84414721,      6.00870609, 
	      5.2993288,      4.69635773,       4.1803813,      3.73364663, 
	      3.3409214,      2.98988771,      2.67116642,      2.37808347, 
	     2.10627484,      1.85321128,      1.61770749,      1.39946342, 
	     1.19866586,      1.01566935,      0.85075897,     0.703995764, 
	    0.575130582,     0.463578522,     0.368436754,     0.288532734, 
	    0.222490266,     0.168802604,     0.125904724,    0.0922389627, 
	   0.0663105547,    0.0467309169,    0.0322484411,    0.0217669662, 
	   0.0143531393,   0.00923430733,   0.00578886038,   0.00353109487, 
	   0.0020927703,   0.00120328926,  0.000670136011,  0.000360891077, 
	 0.000187606522,  9.39672682e-05,   4.5260218e-05,  2.09206464e-05, 
	 9.25993299e-06,  3.91576259e-06,  1.57813565e-06, 
/* Z=84 -> Po */
	              0,      0.90622884,     0.979385972,      1.05837417, 
	     1.14364815,      1.23569608,      1.33504331,      1.44225335, 
	     1.55793178,      1.68272805,      1.81733894,      1.96251142, 
	     2.11904621,      2.28780007,       2.4696908,       2.6657002, 
	     2.87687683,      3.10434175,      3.34929109,      3.61300063, 
	     3.89683008,      4.20222712,      4.53073168,      4.88397932, 
	     5.26370764,      5.67175674,      6.11007643,      6.58072805, 
	     7.08588839,      7.62785101,      8.20903111,      8.83196545, 
	     9.49931335,      10.2138557,      10.9784985,      11.7962608, 
	     12.6702843,      13.6038141,      14.6001997,      15.6628819, 
	     16.7953815,      18.0012798,      19.2842064,      20.6478157, 
	     22.0957527,      23.6316338,      25.2590046,      26.9812984, 
	     28.8017998,      30.7235756,      32.7494354,      34.8818512, 
	     37.1229019,      39.4741745,      41.9366951,      44.5108337, 
	     47.1961937,      49.9915314,      52.8946152,      55.9021378, 
	     59.0095787,      62.2111092,      65.4994583,      68.8657913, 
	     72.2996445,      75.7887878,      79.3191681,      82.8748627, 
	     86.4380264,      89.9889069,      93.5059052,      96.9656219, 
	     100.343071,      103.611809,      106.744293,      109.712173, 
	     112.486771,      115.039589,      117.342941,      119.370628, 
	     121.098778,      122.506638,      123.577599,       124.30011, 
	     124.668663,      124.684837,        124.3582,      123.707138, 
	      122.75956,      121.553329,      120.136543,      118.567299, 
	     116.913261,      115.250549,      113.662323,      112.236649, 
	      111.06398,      110.233879,      109.831497,      109.933441, 
	     110.603432,      111.887939,      113.811836,      116.374451, 
	     119.546349,      123.267113,      127.444359,      131.954529, 
	     136.645477,      141.341095,      145.848282,      149.965759, 
	     153.494965,      156.252258,      158.081985,      158.869583, 
	     158.553665,      157.135925,      154.687988,      151.353683, 
	     147.345963,      142.938019,      138.447861,      134.217194, 
	     130.585297,      127.859558,      126.285156,      126.016579, 
	     127.094231,      129.429413,      132.800323,      136.861145, 
	     141.165222,      145.201233,       148.44043,      150.390427, 
	     150.650375,      148.961105,      145.243607,      139.619904, 
	     132.412476,      124.119637,      115.368446,        106.8489, 
	     99.2366486,      93.1135559,      88.8964005,      86.7838287, 
	     86.7294388,      88.4458923,      91.4400024,      95.0746307, 
	     98.6485748,      101.483269,      103.003983,       102.80426, 
	     100.685631,      96.6683121,      90.9736557,      83.9826965, 
	      76.178566,      68.0818558,      60.1880341,      52.9149208, 
	     46.5656929,      41.3105278,      37.1867409,      34.1151428, 
	     31.9281883,      30.4045982,      29.3049259,      28.4031525, 
	     27.5108261,       26.491806,      25.2674179,      23.8132553, 
	     22.1497822,      20.3293877,      18.4223576,      16.5039864, 
	     14.6441135,      12.8999252,       11.312048,      9.90360069, 
	     8.68151665,      7.63938332,      6.76103401,      6.02427816, 
	      5.4043045,      4.87646437,      4.41831207,      4.01089144, 
	      3.6393435,      3.29296684,      2.96487331,      2.65138793, 
	      2.3513217,      2.06521797,      1.79465115,      1.54163218, 
	     1.30813909,      1.09578598,     0.905621111,     0.738034129, 
	    0.592753649,     0.468908787,     0.365132064,     0.279683918, 
	    0.210582927,     0.155729294,     0.113013938,    0.0804081783, 
	   0.0560319647,    0.0382004417,    0.0254505128,    0.0165496375, 
	   0.0104901353,   0.00647261739,   0.00388203259,    0.0022597427, 
	  0.00127462868,  0.000695509836,  0.000366475899,  0.000186122881, 
	 9.09306764e-05,  4.26456863e-05,  1.91576637e-05,  8.22450693e-06, 
	  3.3660524e-06,  1.30997705e-06,  4.83469023e-07, 
/* Z=85 -> At */
	              0,     0.984225869,      1.06305373,      1.14811301, 
	     1.23988509,      1.33888698,      1.44567406,      1.56084204, 
	      1.6850301,      1.81892383,      1.96325767,      2.11881876, 
	     2.28644991,      2.46705222,      2.66159081,      2.87109566, 
	     3.09666777,      3.33948159,       3.6007905,      3.88192892, 
	     4.18431807,      4.50947046,      4.85899162,      5.23458672, 
	     5.63806391,      6.07133818,      6.53643322,      7.03548861, 
	      7.5707593,      8.14461899,      8.75956345,      9.41821003, 
	     10.1232977,      10.8776913,        11.68437,      12.5464344, 
	     13.4670935,      14.4496632,      15.4975557,      16.6142654, 
	     17.8033638,      19.0684662,      20.4132328,      21.8413219, 
	     23.3563805,      24.9619999,      26.6616821,      28.4587955, 
	     30.3565273,      32.3578262,      34.4653473,      36.6813736, 
	     39.0077515,      41.4458008,      43.9962349,      46.6590576, 
	     49.4334679,      52.3177414,      55.3091393,      58.4037781, 
	     61.5965042,      64.8807983,       68.248642,      71.6904221, 
	     75.1948242,      78.7487488,      82.3372192,      85.9433899, 
	     89.5484695,      93.1317825,      96.6708145,      100.141319, 
	     103.517517,      106.772324,      109.877586,      112.804642, 
	     115.524582,       118.00898,      120.230476,      122.163528, 
	     123.785187,      125.076103,      126.021362,      126.611549, 
	     126.843727,      126.722412,      126.260468,      125.479927, 
	      124.41259,      123.100388,       121.59552,      119.960129, 
	     118.265663,      116.591629,      115.023918,      113.652573, 
	     112.568901,       111.86219,       111.61586,      111.903267, 
	      112.78331,      114.295921,      116.457916,      119.259171, 
	     122.659653,      126.587608,      130.939117,      135.579407, 
	     140.346207,      145.055115,      149.507294,      153.499069, 
	     156.833374,       159.33252,      160.851425,      161.290665, 
	     160.608124,      158.828125,      156.046967,      152.433624, 
	     148.224731,      143.713547,      139.232391,      135.129532, 
	     131.741333,      129.361984,      128.212891,      128.415146, 
	     129.968155,      132.737762,      136.456116,      140.735626, 
	     145.096756,      149.008835,      151.940735,      153.416718, 
	     153.071793,      150.699509,      146.286087,      140.024963, 
	     132.308304,      123.694023,      114.850769,      106.485847, 
	     99.2642212,      93.7286835,       90.231926,      88.8903351, 
	     89.5667648,      91.8857956,      95.2799759,      99.0611572, 
	     102.506828,       104.94912,      105.854012,      104.879921, 
	     101.908569,      97.0453186,      90.5911865,      82.9923935, 
	     74.7761078,      66.4823074,      58.6008606,      51.5218391, 
	     45.5036888,      40.6616821,      36.9755211,      34.3126755, 
	     32.4623833,      31.1741505,      30.1949959,      29.3005619, 
	     28.3169041,      27.1316147,      25.6947727,      24.0116158, 
	     22.1295986,      20.1228638,      18.0767002,      16.0742054, 
	     14.1862497,      12.4653568,      10.9432611,      9.63157463, 
	       8.524683,      7.60394621,      6.84237289,      6.20909309, 
	     5.67315626,      5.20641041,      4.78536463,      4.39211369, 
	      4.0144558,      3.64541721,      3.28238034,      2.92602038, 
	      2.5792017,      2.24596262,      1.93065727,      1.63730383, 
	     1.36914062,      1.12837505,     0.916097224,     0.732320607, 
	    0.576110244,     0.445764333,     0.339019328,     0.253254324, 
	    0.185679793,     0.133498564,    0.0940346569,    0.0648276657, 
	   0.0436940379,    0.0287586078,    0.0184611604,    0.0115431398, 
	  0.00702031422,   0.00414682785,   0.00237534591,   0.00131727161, 
	 0.000706003455,  0.000365032145,  0.000181722426,  8.69283685e-05, 
	 3.98714983e-05,  1.74958623e-05,  7.32744365e-06,  2.92168647e-06, 
	 1.10622625e-06,  3.96636182e-07,  1.34285784e-07, 
/* Z=86 -> Rn */
	              0,      1.06927836,      1.15422595,      1.24583185, 
	     1.34460616,      1.45109558,      1.56588721,      1.68961024, 
	     1.82293975,      1.96659827,       2.1213603,       2.2880547, 
	     2.46756768,      2.66084766,      2.86890697,       3.0928278, 
	      3.3337636,      3.59294581,      3.87168527,      4.17137766, 
	     4.49350786,      4.83965349,       5.2114892,      5.61079025, 
	     6.03943777,      6.49942017,      6.99283981,      7.52191305, 
	       8.088974,      8.69647694,      9.34699821,      10.0432348, 
	     10.7880087,      11.5842609,      12.4350471,      13.3435411, 
	     14.3130226,      15.3468723,      16.4485588,      17.6216278, 
	     18.8696899,      20.1963921,      21.6054096,      23.1004047, 
	     24.6850052,      26.3627701,      28.1371422,      30.0114098, 
	     31.9886475,      34.0716667,      36.2629395,      38.5645485, 
	     40.9780769,      43.5045586,      46.1443634,      48.8971024, 
	     51.7615318,      54.7354431,      57.8155403,      60.9973259, 
	     64.2749939,      67.6413116,      71.0874863,      74.6030884, 
	     78.1759338,      81.7920303,      85.4354858,      89.0884933, 
	     92.7313461,      96.3424301,      99.8983459,      103.374031, 
	      106.74295,      109.977386,      113.048759,      115.928032, 
	     118.586273,      120.995201,      123.127899,      124.959602, 
	     126.468513,      127.636757,       128.45134,      128.905167, 
	     128.998001,      128.737457,      128.139877,      127.231018, 
	     126.046677,      124.632896,      123.045998,      121.352112, 
	     119.626328,      117.951347,      116.415558,       115.11058, 
	      114.12822,      113.556953,      113.477959,      113.960762, 
	     115.058777,      116.804893,      119.207291,      122.245911, 
	     125.869789,      129.995667,      134.508057,      139.261246, 
	     144.083328,      148.782364,        153.1548,      156.995789, 
	     160.111267,      162.331055,      163.522247,      163.602112, 
	     162.549118,      160.411194,      157.309769,      153.438858, 
	     149.057999,      144.478958,      140.046127,      136.111343, 
	     133.004837,      131.003922,      130.302765,       130.98613, 
	     133.010254,      136.194489,      140.225494,      144.675903, 
	     149.037048,      152.763962,      155.329254,      156.280396, 
	     155.294434,      152.222885,      147.120697,      140.253723, 
	     132.082062,      123.218597,      114.366463,       106.24115, 
	     99.4869156,      94.5977859,      91.8542175,       91.284996, 
	     92.6604538,      95.5193329,      99.2260666,      103.050827, 
	     106.261078,      108.211517,      108.419769,      106.617813, 
	     102.772881,      97.0770569,      89.9086685,      81.7732697, 
	     73.2335739,      64.8387833,      57.0626373,      50.2576942, 
	     44.6297913,      40.2341461,      36.9908485,      34.7153664, 
	     33.1579437,       32.045372,      31.1189728,      30.1641884, 
	     29.0290527,      27.6309204,      25.9526958,      24.0311203, 
	     21.9403591,      19.7741222,      17.6290722,      15.5915098, 
	     13.7283268,      12.0824633,      10.6723776,       9.4946909, 
	     8.52888393,      7.74297428,      7.09924841,      6.55932999, 
	      6.0881505,      5.65661097,      5.24296236,      4.83305693, 
	     4.41973686,      4.00164223,      3.58172297,       3.1657002, 
	     2.76063871,      2.37376881,      2.01159739,      1.67932582, 
	     1.38054621,      1.11716759,     0.889516115,     0.696551442, 
	    0.536149144,     0.405405462,     0.300933331,     0.219128102, 
	    0.156389594,     0.109294593,    0.0747195184,     0.049916286, 
	   0.0325474106,    0.0206878241,    0.0128014954,   0.00770090753, 
	   0.0044968687,   0.00254495698,   0.00139357068,  0.000737042748, 
	 0.000375808566,  0.000184375487,  8.68578645e-05,  3.92052971e-05, 
	 1.69169125e-05,   6.9614357e-06,  2.72512216e-06,  1.01213072e-06, 
	 3.55671176e-07,  1.17913743e-07,  3.67676982e-08, 
/* Z=87 -> Fr */
	              0,      1.16211104,      1.25366747,      1.35233819, 
	     1.45866275,      1.57321954,      1.69662869,      1.82955515, 
	     1.97271049,      2.12685728,      2.29281068,      2.47144318, 
	     2.66368723,      2.87053871,      3.09306121,      3.33238912, 
	     3.58973193,      3.86637807,      4.16369915,      4.48315382, 
	     4.82629156,      5.19475794,      5.59029675,      6.01475668, 
	     6.47009134,      6.95836687,      7.48176146,      8.04257107, 
	      8.6432085,      9.28621101,      9.97423267,      10.7100534, 
	     11.4965706,      12.3368034,      13.2338829,      14.1910543, 
	     15.2116613,      16.2991447,      17.4570255,      18.6888962, 
	     19.9983978,      21.3892002,      22.8649864,      24.4294071, 
	     26.0860653,      27.8384686,      29.6899891,      31.6438141, 
	     33.7028923,      35.8698692,      38.1470337,      40.5362167, 
	     43.0387344,      45.6553001,      48.3859138,       51.229763, 
	     54.1851425,      57.2493095,      60.4183922,      63.6872559, 
	      67.049408,      70.4968491,      74.0200119,      77.6076202, 
	     81.2465973,      84.9220428,      88.6171341,      92.3131104, 
	     95.9893188,      99.6232529,      103.190643,      106.665619, 
	     110.020981,      113.228409,      116.258919,       119.08329, 
	     121.672592,      123.998817,       126.03566,      127.759262, 
	     129.149109,      130.188995,      130.868027,      131.181534, 
	     131.132217,      130.730911,      129.997559,      128.961823, 
	     127.663513,       126.15284,      124.490196,      122.745651, 
	     120.997803,      119.332314,      117.839767,      116.612999, 
	     115.743896,      115.319618,      115.418518,      116.105743, 
	     117.428627,      119.412407,      122.056137,      125.329422, 
	     129.170059,      133.483124,      138.141693,      142.989426, 
	     147.845367,      152.510864,      156.778763,      160.444412, 
	     163.318207,       165.23909,      166.087952,      165.800171, 
	     164.376053,      161.887924,      158.482697,      154.379013, 
	     149.858276,      145.248947,      140.904922,      137.178665, 
	     134.390594,      132.797668,      132.563416,      133.733414, 
	     136.218918,      139.792252,      144.095581,      148.664444, 
	     152.965256,      156.444321,      158.584412,      158.962875, 
	     157.304947,      153.524994,      147.749466,      140.316864, 
	     131.752502,      122.718643,       113.94471,       106.14463, 
	     99.9315414,      95.7409973,      93.7735825,       93.966011, 
	     95.9956207,      99.3191223,      103.240471,      106.998848, 
	     109.864082,      111.225891,      110.664558,      107.993668, 
	     103.270302,      96.7731018,      88.9534607,      80.3683472, 
	     71.6054382,      63.2113037,      55.6318626,      49.1720009, 
	     43.9775963,      40.0400887,      37.2202835,      35.2858009, 
	     33.9553032,      32.9423714,      31.9929371,      30.9119873, 
	     29.5775051,       27.941906,      26.0230427,      23.8880653, 
	     21.6338081,      19.3672581,      17.1887016,      15.1793966, 
	     13.3944159,      11.8605051,      10.5781069,      9.52636719, 
	     8.66976643,       7.9651556,      7.36818743,      6.83847141, 
	     6.34309435,       5.8584466,      5.37054205,      4.87414122, 
	     4.37108421,      3.86821628,      3.37525392,      2.90284252, 
	     2.46098089,      2.05789185,      1.69935215,      1.38843095, 
	     1.12557006,     0.908909321,     0.734765351,     0.598181844, 
	    0.493479878,     0.414757639,      0.35630542,      0.31291908, 
	    0.280107528,     0.254202247,     0.232382044,     0.212631524, 
	    0.193652064,     0.174743474,     0.155671984,     0.136537582, 
	    0.117650278,    0.0994221792,    0.0822797343,     0.066598244, 
	   0.0526587293,     0.040626213,    0.0305462852,     0.022356322, 
	   0.0159068704,    0.0109884767,   0.00735960389,   0.00477197999, 
	  0.00299086864,   0.00180900819,   0.00105409417, 
/* Z=88 -> Ra */
	              0,       1.2634877,      1.36218297,      1.46847892, 
	     1.58294618,      1.70619655,      1.83888495,      1.98171318, 
	       2.135432,      2.30084443,       2.4788096,      2.67024469, 
	      2.8761301,      3.09751153,      3.33550525,      3.59130025, 
	     3.86616302,      4.16144228,      4.47857141,       4.8190732, 
	     5.18456459,      5.57676029,       5.9974761,      6.44863224, 
	     6.93225956,      7.45050001,      8.00561047,      8.59996605, 
	     9.23606014,      9.91650963,      10.6440487,      11.4215374, 
	     12.2519484,      13.1383762,      14.0840263,       15.092206, 
	     16.1663265,      17.3098793,      18.5264359,       19.819622, 
	     21.1931057,      22.6505718,      24.1956921,      25.8321056, 
	     27.5633698,       29.392931,        31.32407,      33.3598595, 
	     35.5030975,      37.7562523,      40.1213875,      42.6000824, 
	     45.1933517,      47.9015503,      50.7242889,      53.6603127, 
	     56.7074013,      59.8622665,      63.1204147,       66.476059, 
	     69.9219666,       73.449379,       77.047905,      80.7053909, 
	     84.4078674,      88.1395035,      91.8825226,      95.6172333, 
	     99.3220444,      102.973549,      106.546638,      110.014702, 
	     113.349869,      116.523354,      119.505821,      122.267937, 
	     124.780891,      127.017113,      128.950989,      130.559753, 
	     131.824326,      132.730377,      133.269196,      133.438812, 
	     133.244934,      132.701859,      131.833206,      130.672607, 
	     129.264023,      127.661812,       125.93042,      124.143707, 
	     122.383636,      120.738602,      119.301048,       118.16465, 
	      117.42086,      117.155029,      117.442131,        118.3423, 
	     119.896233,      122.120949,      125.005882,      128.509872, 
	     132.559265,       137.04747,      141.836212,      146.758957, 
	     151.626343,      156.234009,      160.372345,      163.838287, 
	     166.448257,      168.051895,      168.545456,      167.883865, 
	     166.090363,      163.262344,      159.572235,      155.262833, 
	     150.636017,      146.034958,      141.820404,      138.342056, 
	      135.90712,      134.748535,      134.996033,      136.653442, 
	     139.585648,      143.517883,      148.049347,      152.681656, 
	     156.860977,      160.030807,      161.690598,      161.454086, 
	     159.100494,      154.611191,      148.186111,      140.235611, 
	     131.346466,      122.223808,      113.614601,        106.2211, 
	     100.615112,      97.1646729,      95.9839783,      96.9146652, 
	     99.5422134,      103.246567,      107.280159,       110.86264, 
	     113.278702,      113.965103,      112.574631,      109.009041, 
	      103.41761,      96.1636734,      87.7656784,      78.8227921, 
	     69.9364319,      61.6387215,      54.3367424,       48.278965, 
	     43.5458183,      40.0634079,      37.6356354,      35.9882278, 
	     34.8167648,        33.83144,      32.7925491,      31.5329876, 
	     29.9666557,      28.0839691,      25.9374599,      23.6214619, 
	     21.2499447,      18.9361362,      16.7763863,      14.8398218, 
	      13.163909,      11.7553682,       10.595109,      9.64569187, 
	     8.85966682,      8.18750095,      7.58407307,      7.01318979, 
	     6.44994068,       5.8810792,      5.30381823,      4.72356796, 
	     4.15114164,      3.59991813,      3.08331442,      2.61282229, 
	     2.19670725,       1.8393805,      1.54136729,       1.2997371, 
	     1.10884643,     0.961239219,     0.848568201,     0.762429595, 
	    0.695038617,     0.639707386,     0.591119826,     0.545419097, 
	    0.500143647,      0.45405212,     0.406881899,     0.359078974, 
	    0.311532408,     0.265336424,     0.221594453,     0.181273445, 
	    0.145108804,     0.113556817,    0.0867876485,    0.0647103414, 
	   0.0470197536,    0.0332558155,    0.0228661299,    0.0152643546, 
	  0.00987898558,   0.00618935144,   0.00374791608,   0.00218988117, 
	  0.00123245502,  0.000666857988,  0.000346223096, 
/* Z=89 -> Ac */
	              0,      1.37422657,      1.48063254,      1.59515774, 
	      1.7184062,      1.85102391,      1.99370289,      2.14718318, 
	      2.3122561,      2.48976779,       2.6806221,      2.88578391, 
	     3.10628271,      3.34321618,      3.59775496,      3.87114453, 
	     4.16471052,      4.47986317,      4.81809998,      5.18100977, 
	     5.57027817,      5.98769045,      6.43513441,      6.91460657, 
	     7.42821312,      7.97817326,      8.56682301,      9.19661808, 
	     9.87013054,      10.5900574,      11.3592119,      12.1805305, 
	      13.057065,      13.9919777,      14.9885454,      16.0501385, 
	     17.1802216,      18.3823395,      19.6601028,      21.0171642, 
	     22.4572124,      23.9839287,      25.6009789,      27.3119583, 
	     29.1203747,      31.0295925,      33.0427895,      35.1629028, 
	     37.3925629,      39.7340279,      42.1891251,      44.7591476, 
	     47.4447823,      50.2460098,      53.1620102,      56.1910515, 
	     59.3303833,      62.5761299,      65.9231567,       69.364975, 
	     72.8936005,      76.4994888,      80.1713791,      83.8962402, 
	     87.6591949,      91.4434509,      95.2303009,      98.9990921, 
	     102.727333,      106.390732,      109.963379,      113.417953, 
	     116.725998,      119.858292,      122.785271,      125.477577, 
	     127.906662,      130.045486,      131.869308,      133.356598, 
	     134.489883,      135.256821,      135.651154,      135.673721, 
	     135.333435,      134.648132,      133.645294,      132.362564, 
	     130.848099,        129.1604,      127.367912,      125.548141, 
	     123.786232,      122.172989,      120.802406,      119.768539, 
	     119.161942,       119.06559,      119.550499,      120.671188, 
	     122.461159,      124.928642,      128.053024,      131.782059, 
	     136.030487,      140.680008,      145.581375,      150.558304, 
	      155.41391,      159.939072,       163.92308,      167.165894, 
	     169.491592,      170.762131,       170.89064,      169.852905, 
	     167.696045,      164.542892,      160.591293,      156.107269, 
	     151.411545,      146.859604,      142.816055,      139.624405, 
	     137.574966,      136.873138,      137.611786,      139.750916, 
	     143.108017,      147.361862,      152.070908,      156.706757, 
	     160.700653,      163.499573,      164.626434,      163.737854, 
	     160.672318,      155.481674,      148.440247,      140.028305, 
	     130.889114,      121.763046,      113.405098,      106.495422, 
	     101.554489,      98.8746948,      98.4786224,      100.111504, 
	      103.26992,      107.264038,      111.304932,      114.604752, 
	     116.475456,      116.411797,      114.147324,      109.676361, 
	     103.240593,      95.2840424,       86.384697,       77.173996, 
	     68.2560959,      60.1379585,      53.1788864,      47.5647812, 
	     43.3079224,      40.2699852,      38.2021942,      36.7949677, 
	     35.7287064,      34.7181473,      33.5447121,      32.0739517, 
	     30.2579231,      28.1247406,      25.7591381,      23.2786503, 
	     20.8096848,      18.4671688,      16.3398972,      14.4826536, 
	     12.9148092,      11.6242771,      10.5751228,      9.71698475, 
	     8.99453545,      8.35564232,      7.75733614,      7.16921568, 
	     6.57437515,      5.96828794,      5.35624361,      4.75005531, 
	     4.16462994,      3.61491132,      3.11352015,      2.66924095, 
	     2.28636265,      1.96475232,      1.70049012,      1.48684847, 
	     1.31541324,      1.17717147,      1.06344068,     0.966560423, 
	    0.880323827,     0.800157785,     0.723094165,     0.647587061, 
	    0.573233783,     0.500454366,     0.430174053,     0.363538861, 
	    0.301685333,     0.245572329,     0.195875674,     0.152939707, 
	     0.11677672,    0.0871025845,    0.0633964837,    0.0449731275, 
	   0.0310570225,    0.0208506193,    0.0135903647,   0.00858742371, 
	  0.00525224721,   0.00310436171,   0.00177011138,  0.000971955771, 
	  0.00051296223,  0.000259686552,  0.000125842023, 
/* Z=90 -> Th */
	              0,      1.49529815,      1.61003745,      1.73344958, 
	     1.86617279,      2.00888968,      2.16233039,      2.32727575, 
	     2.50455928,      2.69507194,       2.8997643,      3.11965013, 
	     3.35581088,      3.60939789,      3.88163757,      4.17383432, 
	     4.48737526,      4.82373285,      5.18447018,      5.57124519, 
	     5.98581314,      6.43003082,      6.90586042,      7.41537428, 
	     7.96075535,      8.54430199,      9.16842937,      9.83567047, 
	     10.5486813,      11.3102341,      12.1232214,      12.9906549, 
	      13.915658,      14.9014664,       15.951416,      17.0689392, 
	     18.2575531,      19.5208435,      20.8624516,      22.2860565, 
	     23.7953453,      25.3939991,      27.0856495,      28.8738499, 
	     30.7620335,      32.7534714,      34.8512154,      37.0580444, 
	     39.3764038,      41.8083267,      44.3553658,      47.0185051, 
	     49.7980804,      52.6936646,      55.7039795,      58.8267899, 
	     62.0587769,      65.3954468,      68.8309937,      72.3582001, 
	     75.9683075,      79.6509323,      83.3939514,      87.1834259, 
	     91.0035477,      94.8365707,      98.6628342,       102.46077, 
	     106.206947,      109.876259,      113.442009,      116.876228, 
	      120.14994,      123.233566,      126.097412,      128.712204, 
	     131.049744,       133.08371,      134.790375,      136.149582, 
	     137.145645,      137.768417,      138.014236,      137.886917, 
	     137.398727,      136.571167,      135.435699,      134.034103, 
	     132.418671,      130.652054,      128.806625,      126.963371, 
	     125.210373,      123.640572,      122.349091,      121.429939, 
	     120.972267,      121.056099,      121.747917,      123.096039, 
	     125.126137,      127.837219,      131.198135,      135.145355, 
	     139.581802,       144.37767,      149.373001,      154.382446, 
	     159.202347,      163.620056,      167.425064,      170.421906, 
	     172.443863,      173.366837,      173.122314,      171.708115, 
	      169.19603,       165.73465,      161.546829,      156.920715, 
	     152.194046,      147.732117,      143.900177,      141.032104, 
	     139.397675,      139.171387,      140.406403,      143.017075, 
	     146.773071,      151.307678,      156.141357,      160.719971, 
	     164.465485,      166.834793,      167.380966,      165.809616, 
	     162.023331,      156.147156,      148.529922,      139.718872, 
	      130.40802,      121.364738,      113.342133,        106.9879, 
	     102.761627,       100.87262,      101.247772,      103.535912, 
	     107.148972,      111.335419,      115.275948,      118.187729, 
	     119.422218,      118.542709,      115.370811,      109.996758, 
	     102.753532,      94.1605377,      84.8466873,      75.4648056, 
	     66.6101074,      58.7534676,      52.1976013,       47.060154, 
	     43.2834015,      40.6665764,      38.9135551,       37.687233, 
	     36.6617584,      35.5651474,      34.2073021,      32.4915314, 
	     30.4104042,      28.0293503,      25.4625835,       22.846447, 
	     20.3145943,      17.9784908,      15.9149876,      14.1614923, 
	     12.7178574,      11.5534229,      10.6170483,      9.84803486, 
	     9.18615246,      8.57944775,      7.98918104,      7.39176989, 
	      6.7781477,      6.15120697,      5.52214813,      4.90653181, 
	     4.32065678,      3.77872562,      3.29099298,      2.86292315, 
	     2.49522948,      2.18456268,      1.92459655,      1.70725715, 
	     1.52389836,      1.36627245,      1.22722602,      1.10109913, 
	    0.983856142,     0.873001933,     0.767350674,     0.666714966, 
	    0.571571648,     0.482747942,     0.401156694,     0.327595264, 
	    0.262612402,     0.206437767,     0.158966437,     0.119785048, 
	   0.0882275179,    0.0634471476,    0.0444943011,    0.0303900205, 
	   0.0201888755,    0.0130267879,   0.00815197825,    0.0049398304, 
	  0.00289382041,   0.00163603225,  0.000891013944,  0.000466574129, 
	 0.000234438543,   0.00011279602,  5.18505076e-05, 
/* Z=91 -> Pa */
	              0,      1.62769306,      1.75143683,      1.88444364, 
	     2.02738786,      2.18098998,       2.3460207,      2.52330232, 
	     2.71371365,      2.91819119,      3.13773489,      3.37340951, 
	     3.62634921,      3.89776158,      4.18893051,      4.50122166, 
	     4.83608389,      5.19505548,      5.57976723,      5.99194479, 
	     6.43341589,      6.90611029,      7.41206646,      7.95343208, 
	     8.53246784,      9.15155125,      9.81317711,      10.5199566, 
	     11.2746248,      12.0800323,       12.939146,       13.855052, 
	     14.8309441,      15.8701229,      16.9759846,      18.1520138, 
	     19.4017715,      20.7288818,      22.1370106,      23.6298447, 
	     25.2110691,      26.8843422,      28.6532574,      30.5213051, 
	     32.4918289,      34.5679893,      36.7526932,       39.048542, 
	     41.4577675,       43.982151,      46.6229553,      49.3808289, 
	     52.2557106,      55.2467499,      58.3521805,      61.5692215, 
	     64.8939667,      68.3212662,      71.8446121,      75.4560318, 
	     79.1459579,      82.9031601,      86.7146225,      90.5654755, 
	     94.4389725,      98.3164062,      102.177177,      105.998787, 
	     109.756943,      113.425682,      116.977623,      120.384178, 
	     123.615929,      126.643028,      129.435745,      131.964996, 
	      134.20314,      136.124619,      137.706909,      138.931396, 
	     139.784363,      140.258011,       140.35141,      140.071548, 
	     139.434189,      138.464645,       137.19838,      135.681335, 
	     133.970047,       132.13118,      130.240891,      128.383514, 
	     126.649757,      125.134369,      123.933189,      123.139717, 
	     122.841133,      123.113937,      124.019478,      125.599297, 
	     127.870728,      130.823013,      134.414215,      138.569183, 
	     143.179153,      148.102966,      153.170395,      158.187775, 
	     162.945694,      167.229156,      170.829254,      173.556625, 
	     175.255386,      175.816864,      175.192169,      173.402222, 
	     170.543915,       166.79155,      162.392273,      157.654877, 
	     152.931976,      148.595901,      145.009415,      142.493378, 
	     141.293777,      141.551483,      143.278137,      146.341797, 
	     150.465149,      155.238556,      160.148575,      164.620804, 
	     168.074326,      169.982666,      169.935181,      167.691284, 
	     163.220276,      156.720184,       148.61087,      139.499939, 
	     130.122757,      121.262428,      113.658379,      107.914948, 
	     104.422531,      103.302902,      104.387703,      107.234848, 
	     111.181686,      115.428497,      119.140617,      121.554649, 
	     122.073372,      120.335907,      116.253891,      110.009689, 
	     102.018791,      92.8637238,      83.2104187,      73.7200165, 
	     64.9685364,      57.3852768,      51.2165108,      46.5177231, 
	     43.1723518,      40.9317436,      39.4680748,      38.4308929, 
	     37.4983826,        36.41642,      35.0213165,      33.2453957, 
	     31.1074562,      28.6922817,       26.124279,      23.5405598, 
	     21.0674973,      18.8037853,      16.8111172,       15.112381, 
	     13.6959848,      12.5244141,      11.5447502,      10.6991816, 
	     9.93392754,      9.20558929,      8.48459816,      7.75590181, 
	     7.01744652,      6.27718401,      5.54938364,      4.85093784, 
	     4.19815779,       3.6043992,      3.07861257,       2.6248033, 
	     2.24224448,      1.92623413,      1.66917205,      1.46174908, 
	     1.29407275,      1.15661395,      1.04090607,     0.939977705, 
	    0.848536491,     0.762946367,       0.6810534,     0.601917028, 
	    0.525498211,     0.452347428,     0.383321434,     0.319349587, 
	    0.261257708,     0.209650651,     0.164850265,     0.126880035, 
	   0.0954864547,    0.0701863766,    0.0503294282,    0.0351656675, 
	   0.0239102636,    0.0157989822,    0.0101305135,   0.00629413733, 
	  0.00378314848,   0.00219612126,   0.00122908072,  0.000661943981, 
	 0.000342399173,  0.000169755906,  8.04938245e-05, 
/* Z=92 -> U */
	              0,      1.77268696,      1.90617299,      2.04955173, 
	     2.20353532,       2.3688848,      2.54641223,      2.73698449, 
	     2.94152594,      3.16102242,      3.39652371,      3.64914799, 
	     3.92008471,      4.21059942,      4.52203512,      4.85582018, 
	     5.21346855,      5.59658432,      6.00686693,      6.44611406, 
	     6.91622591,      7.41920757,      7.95717192,      8.53234482, 
	     9.14706612,      9.80379009,      10.5050917,      11.2536631, 
	     12.0523119,      12.9039679,      13.8116732,      14.7785816, 
	     15.8079548,      16.9031544,      18.0676327,      19.3049183, 
	     20.6186161,      22.0123711,      23.4898682,      25.0547924, 
	     26.7108173,       28.461565,      30.3105717,      32.2612495, 
	     34.3168411,      36.4803696,      38.7545815,      41.1418762, 
	     43.6442451,      46.2632027,      48.9996834,      51.8539696, 
	     54.8255959,      57.9132309,      61.1145973,      64.4263458, 
	     67.8439407,      71.3615494,      74.9719315,      78.6663208, 
	     82.4343338,      86.2638397,      90.1409302,      94.0497971, 
	     97.9727249,      101.890068,      105.780281,      109.619942, 
	     113.383911,      117.045448,        120.5765,      123.947952, 
	     127.129997,      130.092636,      132.806183,      135.241913, 
	     137.372818,      139.174316,      140.625198,      141.708603, 
	     142.412933,       142.73291,      142.670563,      142.236145, 
	     141.449097,      140.338638,      138.944336,      137.316376, 
	     135.515427,      133.612183,      131.686447,      129.825653, 
	     128.122894,      126.674339,      125.576149,      124.920776, 
	       124.7929,      125.264938,      126.392502,      128.209732, 
	     130.725082,      133.917618,      137.734207,      142.087982, 
	     146.858414,       151.89328,      157.012634,      162.015167, 
	     166.686783,      170.811371,      174.183136,      176.620361, 
	     177.979584,      178.169144,      177.161133,      175.000458, 
	     171.809738,       167.78891,      163.208588,      158.396744, 
	     153.718643,      149.550797,      146.250137,      144.120895, 
	     143.381729,      144.136749,      146.353729,      149.853195, 
	      154.31102,      159.276398,      164.205307,      168.507843, 
	     171.606064,      172.996506,      172.311066,      169.368347, 
	     164.208115,      157.102936,      148.542648,       139.19133, 
	     129.819351,      121.217216,      114.101334,       109.02346, 
	     106.296913,      105.950577,      107.718811,      111.070343, 
	     115.273422,      119.488617,       122.87648,      124.704224, 
	     124.436371,      121.796341,      116.791107,      109.696899, 
	     101.009773,       91.369957,      81.4724121,      71.9768524, 
	     63.4299278,      56.2097282,      50.4981804,      46.2830772, 
	     43.3862801,      41.5113487,      40.3013268,      39.3967857, 
	     38.4854126,      37.3366966,      35.8187637,      33.8976593, 
	     31.6222038,      29.0994091,      26.4660358,      23.8614826, 
	     21.4058914,      19.1858215,      17.2480412,      15.6007233, 
	     14.2201681,      13.0608997,       12.066762,      11.1811361, 
	     10.3549252,      9.55161953,      8.74935913,      7.94041157, 
	     7.12874508,      6.32651854,      5.55024147,      4.81723022, 
	     4.14276886,      3.53819299,      3.00992084,      2.55933094, 
	     2.18329406,       1.8751384,      1.62582195,      1.42512417, 
	     1.26271009,      1.12897515,      1.01562965,     0.916021585, 
	    0.825228989,     0.739968657,     0.658379018,     0.579729676, 
	    0.504105926,     0.432104319,     0.364565194,     0.302356541, 
	    0.246216059,     0.196649894,     0.153883606,     0.117856257, 
	   0.0882478654,    0.0645296127,    0.0460263453,    0.0319823176, 
	   0.0216224026,    0.0142034264,   0.00905210339,   0.00558871683, 
	  0.00333721889,   0.00192413072,   0.00106927846,  0.000571666285, 
	  0.00029345206,  0.000144337057,  6.78769429e-05, 
/* Z=93 -> Np */
	              0,       1.9315691,      2.07559729,       2.2301898, 
	     2.39609885,       2.5741272,      2.76513076,      2.97002244, 
	     3.18977499,       3.4254241,      3.67807269,      3.94889379, 
	     4.23913383,      4.55011797,      4.88325167,      5.24002647, 
	     5.62202263,      6.03091288,      6.46846771,      6.93655682, 
	     7.43715382,      7.97233963,      8.54430389,      9.15535164, 
	     9.80789948,      10.5044813,      11.2477503,      12.0404749, 
	     12.8855419,      13.7859516,      14.7448206,      15.7653685, 
	     16.8509197,      18.0048885,      19.2307816,      20.5321674, 
	     21.9126797,      23.3759823,      24.9257622,      26.5657005, 
	     28.2994366,      30.1305447,      32.0624924,      34.0985947, 
	     36.2419739,      38.4954987,       40.861721,      43.3428268, 
	     45.9405479,       48.656086,      51.4900284,      54.4422684, 
	     57.5118904,      60.6970673,      63.9949722,      67.4016495, 
	      70.911911,      74.5192032,      78.2155228,      81.9912949, 
	     85.8352737,      89.7344284,      93.6739197,       97.637001, 
	      101.60498,      105.557266,      109.471382,      113.323021, 
	     117.086197,      120.733467,      124.236153,      127.564659, 
	     130.688904,      133.578827,      136.204926,      138.538956, 
	     140.554626,      142.228546,      143.541016,      144.477036, 
	     145.027344,      145.189362,      144.968185,      144.377594, 
	     143.440735,      142.190872,      140.671768,      138.937836, 
	     137.053864,      135.094467,      133.142944,      131.289566, 
	     129.629517,      128.260025,      127.277046,      126.771523, 
	     126.825043,      127.505348,      128.861679,      130.920258, 
	     133.680145,      137.109772,      141.144547,      145.685776, 
	     150.601334,      155.728195,      160.877304,      165.840637, 
	     170.400513,      174.341034,      177.460999,      179.588013, 
	     180.592575,      180.401505,      179.009033,       176.48494, 
	     172.977798,      168.712814,      163.983017,      159.133789, 
	     154.540573,      150.580994,      147.602707,      145.889786, 
	     145.630325,      146.888931,      149.587738,      153.499298, 
	     158.253723,       163.36171,      168.252823,      172.327209, 
	     175.016083,      175.845474,      174.495682,      170.848984, 
	     165.018082,      157.349823,      148.400955,      138.886154, 
	     129.602463,      121.338074,      114.776237,      110.407059, 
	      108.45993,      108.867409,      111.267159,       115.04361, 
	     119.404358,      123.481468,      126.442993,      127.598618, 
	     126.484421,      122.914459,       116.99337,      109.089188, 
	     99.7725677,      89.7325592,      79.6825409,      70.2700348, 
	     62.0028381,      55.2008972,      49.9781494,      46.2546959, 
	     43.7941437,      42.2581253,       41.267746,      40.4618645, 
	     39.5437088,      38.3103561,       36.663208,      34.6009979, 
	     32.1995621,      29.5839691,        26.89884,      24.2819824, 
	     21.8446369,      19.6601353,      17.7608128,      16.1419411, 
	     14.7704382,       13.595911,      12.5618248,      11.6150246, 
	     10.7125359,      9.82526016,       8.9387598,      8.05171394, 
	      7.1728301,      6.31707335,      5.50187302,      4.74388218, 
	     4.05657291,      3.44879079,      2.92423201,      2.48167348, 
	     2.11575174,      1.81806386,      1.57837808,      1.38578916, 
	     1.22969747,      1.10054636,     0.990296066,     0.892647803, 
	    0.803057194,     0.718587637,      0.63765955,     0.559745431, 
	    0.485053837,     0.414233953,      0.34812206,     0.287542045, 
	    0.233162999,     0.185412169,     0.144437149,      0.11010842, 
	   0.0820520893,    0.0597027577,    0.0423662104,    0.0292835627, 
	   0.0196896046,    0.0128605133,   0.00814807229,   0.00499988813, 
	  0.00296669593,   0.00169924006,   0.00093783997,  0.000497826666, 
	  0.00025365522,  0.000123800535,  5.77517021e-05, 
/* Z=94 -> Pu */
	              0,      2.10576892,      2.26120639,      2.42792511, 
	     2.60671854,      2.79843283,      3.00397015,      3.22429109, 
	     3.46041846,      3.71344066,      3.98451495,      4.27487135, 
	     4.58581543,      4.91873312,      5.27509451,      5.65645647, 
	     6.06446695,      6.50086927,      6.96750546,      7.46631813, 
	     7.99935627,      8.56877804,      9.17684937,      9.82595253, 
	     10.5185852,      11.2573566,      12.0449991,      12.8843575, 
	     13.7783909,      14.7301741,      15.7428865,      16.8198166, 
	     17.9643421,      19.1799316,      20.4701309,      21.8385448, 
	     23.2888241,      24.8246441,      26.4496899,      28.1676083, 
	     29.9820061,      31.8963909,      33.9141426,      36.0384598, 
	     38.2723236,      40.6184273,      43.0791168,      45.6563225, 
	       48.35149,       51.165493,      54.0985451,      57.1501007, 
	     60.3187752,      63.6022186,      66.9970093,      70.4985657, 
	     74.1009903,      77.7970123,      81.5778275,      85.4330063, 
	     89.3504257,      93.3161392,      97.3143692,      101.327393, 
	     105.335579,      109.317368,      113.249382,       117.10643, 
	     120.861771,      124.487259,      127.953659,      131.231003, 
	     134.289017,      137.097672,      139.627777,      141.851654, 
	     143.744034,      145.282745,      146.449783,      147.232239, 
	     147.623306,      147.623291,       147.24057,      146.492477, 
	     145.406097,      144.018784,       142.37854,      140.543976, 
	     138.583984,      136.576965,      134.609467,      132.774399, 
	     131.168625,      129.890015,       129.03389,      128.689087, 
	     128.933533,      129.829712,      131.419891,      133.721832, 
	     136.724686,      140.385956,      144.629333,       149.34433, 
	     154.387344,      159.585083,      164.740067,      169.638535, 
	     174.060425,      177.791473,      180.636642,      182.434387, 
	      183.07077,      182.492371,       180.71669,      177.838867, 
	     174.033508,        169.5504,      164.703644,      159.853943, 
	     155.384216,      151.669983,      149.046295,      147.773666, 
	     148.006638,      149.768173,      152.933777,      157.228348, 
	     162.238266,      167.439224,      172.239487,      176.035126, 
	     178.273056,      178.514954,      176.494644,      172.161346, 
	      165.70137,      157.533722,      148.276764,      138.687851, 
	     129.580887,      121.731018,      115.778152,      112.142487, 
	     110.964539,      112.079781,      115.032951,      119.132179, 
	     123.535812,      127.360527,      129.795013,      130.202713, 
	     128.198914,      123.691643,      116.881554,      108.222397, 
	     98.3499222,      87.9906464,      77.8649597,      68.5985489, 
	     60.6539803,      54.2908173,      49.5572395,      46.3117867, 
	     44.2686653,      43.0574608,      42.2864227,      41.5992813, 
	     40.7177048,      39.4650078,      37.7705383,       35.657341, 
	     33.2182617,      30.5864201,      27.9058819,      25.3071861, 
	     22.8904305,      20.7169609,      18.8089237,      17.1550331, 
	     15.7200594,      14.4557352,      13.3110132,      12.2402496, 
	     11.2085924,      10.1944351,      9.18933105,       8.1960144, 
	     7.22528887,      6.29253626,      5.41439295,      4.60602522, 
	      3.8791914,      3.24114466,      2.69431019,       2.2365725, 
	     1.86200345,      1.56184304,      1.32556605,      1.14190567, 
	    0.999739528,      0.88878566,     0.800084352,      0.72627157, 
	     0.66166836,     0.602220535,     0.545329809,     0.489615589, 
	    0.434643775,     0.380652279,      0.32829529,     0.278422028, 
	    0.231898978,     0.189478606,     0.151714534,     0.118918508, 
	   0.0911534205,    0.0682543591,    0.0498694293,    0.0355118513, 
	   0.0246154536,    0.0165870339,    0.0108506307,   0.00688069919, 
	  0.00422311155,    0.0025046668,     0.001432986,  0.000789453217, 
	 0.000418003474,  0.000212293031,  0.000103200618, 
/* Z=95 -> Am */
	              0,      2.29702711,      2.46482611,      2.64467096, 
	     2.83739877,      3.04390144,      3.26512885,      3.50209165, 
	     3.75586462,      4.02759027,      4.31848335,      4.62983179, 
	     4.96300268,      5.31944561,      5.70069551,      6.10837698, 
	     6.54420757,      7.01000214,      7.50767469,       8.0392437, 
	     8.60683441,      9.21268082,      9.85912704,      10.5486336, 
	     11.2837763,      12.0672417,      12.9018373,      13.7904835, 
	     14.7362099,      15.7421608,      16.8115807,       17.947813, 
	     19.1542912,      20.4345284,      21.7921047,      23.2306461, 
	     24.7538204,      26.3652992,      28.0687447,       29.867775, 
	     31.7659321,      33.7666473,      35.8731956,      38.0886421, 
	     40.4157982,      42.8571663,      45.4148483,      48.0905113, 
	     50.8852806,      53.7996674,       56.833477,       59.985714, 
	     63.2544746,      66.6368484,      70.1288071,        73.72509, 
	     77.4190826,      81.2027359,      85.0664139,      88.9988403, 
	     92.9869766,      97.0159454,      101.068993,      105.127441, 
	     109.170708,        113.1763,      117.119911,      120.975548, 
	     124.715714,      128.311646,      131.733627,       134.95137, 
	     137.934525,      140.653229,       143.07872,      145.184082, 
	     146.945084,      148.341095,      149.355911,      149.978897, 
	     150.205902,      150.040268,      149.493805,      148.587616, 
	     147.352768,      145.830826,      144.074036,      142.145203, 
	     140.117279,       138.07225,      136.099747,      134.294983, 
	     132.756195,      131.581406,      130.864807,      130.692566, 
	     131.138397,       132.25885,      134.088684,      136.636627, 
	     139.881485,      143.769424,      148.212326,       153.08783, 
	     158.241196,      163.489304,      168.627075,      173.435974, 
	     177.694839,       181.19252,      183.741684,      185.193314, 
	     185.450562,      184.481125,      182.326752,      179.108551, 
	     175.027039,      170.355896,      165.428848,      160.619705, 
	      156.31604,      152.888016,      150.654373,      149.848618, 
	     150.588516,      152.853027,      156.469604,      161.115677, 
	     166.335373,      171.572723,      176.219147,      179.672485, 
	     181.401642,      181.010147,      178.290817,      173.263626, 
	     166.190323,      157.561127,       148.05246,      138.458252, 
	     129.601166,      122.234154,      116.944679,      114.075233, 
	     113.671997,      115.471184,      118.926224,      123.274178, 
	     127.632469,      131.112762,      132.935471,      132.528107, 
	     129.593613,      124.139297,      116.463493,        107.1036, 
	     96.7549362,      86.1732025,      76.0758896,      67.0563431, 
	     59.5220642,       53.664135,       49.459362,      46.7020187, 
	     45.0571365,       44.125042,      43.5058517,      42.8540306, 
	     41.9159622,      40.5472755,      38.7106743,      36.4580879, 
	     33.9029999,      31.1893234,      28.4624271,      25.8466358, 
	     23.4311066,      21.2645988,      19.3577003,      17.6905346, 
	     16.2232456,      14.9070196,      13.6937294,      12.5430536, 
	     11.4266481,      10.3294353,      9.24860859,       8.1910696, 
	     7.17006493,      6.20173788,      5.30205393,      4.48444653, 
	     3.75828862,       3.1281836,      2.59395981,      2.15119195, 
	     1.79207146,      1.50644743,      1.28288746,      1.10964942, 
	    0.975488007,     0.870258689,     0.785309434,     0.713674545, 
	    0.650098801,      0.59092921,     0.533914685,     0.477950931, 
	    0.422803551,      0.36883691,     0.316767305,     0.267453998, 
	    0.221734866,      0.18030861,     0.143661946,     0.112036899, 
	   0.0854317695,    0.0636277497,    0.0462329276,    0.0327353962, 
	   0.0225580167,    0.0151086878,   0.00982174743,   0.00618795538, 
	  0.00377248321,   0.00222186767,   0.00126203301,  0.000690074987, 
	 0.000362547347,  0.000182642703,  8.80420339e-05};
                          
numero *Rho_density(int Z)  {return Rho_density_data+(Z-1)*Rho_imax;}
#ifdef jga_spline
void Rho_density(spline &rr, int Z)
  {rr.init(Rho_radius,Rho_density(Z),Rho_imax);}
numero Rho_rmt(int Z, numero valmin)
{
  spline density;  Rho_density(density, Z);
  numero val=0, r=1.5, rmax=Rho_radius[Rho_imax-1];
  int    Zt;  if(Z<10) Zt=Z; else Zt=10;
  while(val<valmin && r<rmax) {r+=0.1; val=(density.integ(0,r)-(Z-Zt))/Zt;}
  density.free();
  return r;
}
numero Rho_rmt(int Z) {return Rho_rmt(Z,0.8);}
#ifdef jga_Poisson
#ifdef jga_Salvat
#ifdef jga_exchange
void Rho_to_potential(int Z, spline &rho, spline &rho_1,
                      numero *r, numero *rV)
{
  numero val,val_1;
  int i;
  for(i=0; i<rho.n; i++) {       
    r[i]=rho.x[i];               
    rV[i]=Poisson(r[i],rho_1)-Z;
    if(i>0) {
      val=rho.y[i]/(4*pi*sqr(r[i]));
      val_1=rho_1.y[i]/(4*pi*sqr(r[i]));
      rV[i] += r[i] * Vxc_Barth_Hedin(val_1,val_1-val/2);
  } }
  rV[rho.n-1]=0;
}
numero Rho_dbound(int Z, numero *r, numero *P, numero *Q,
                  numero &E, int n, int k, int &nmax)
{
  Salvat state;          
  spline rho,rho_1;      
  numero rV[Rho_imax];   
  int i, it, nit=9;      
  numero E0;
  Rho_density(rho, Z);                
  Rho_to_potential(Z,rho,rho,r,rV);   
  for(it=0; it<nit; it++) {
    if(it>0) {
      rho_1.init(rho.x,rho.y,rho.n);
      for(i=0; i<rho.n; i++)
        rho_1.y[i]=rho.y[i]-(sqr(state.P[i])+sqr(state.Q[i]));
      rho_1.init();
      Rho_to_potential(Z,rho,rho_1,r,rV);  rho_1.free();
    }
    state.vint(r,rV,rho.n);
    E0=E;  state.dbound(E, 0.000001,n,k);
    
  }
  nmax=state.nmax;
  for(i=0; i<state.nmax; i++) {
    r[i]=state.rr[i];
    P[i]=state.P[i];  Q[i]=state.Q[i];
  }
  rho.free();  state.free();
  return (E-E0)/E;
}
numero Rho_dbound(int Z, spline &P, spline &Q, numero &E, int n, int k)
{
  numero *r, *PP, *QQ, err;
  int    i, nmax;
  r=new numero [Rho_imax];
  PP=new numero [Rho_imax];  QQ=new numero [Rho_imax];
  err=Rho_dbound(Z,r,PP,QQ,E,n,k,nmax);
  P.free();  Q.free();  P.alloc(nmax);  Q.alloc(nmax);
  for(i=0; i<nmax; i++)  {P.put(i,r[i],PP[i]);  Q.put(i,r[i],QQ[i]);}
  P.init();  Q.init();
  delete [] r;  delete [] PP;  delete [] QQ;
  return err;
}
void Rho_dfree(int Z, numero *r, numero *P, numero *Q,
               numero E, int k, int &nmax)
{
  Salvat state;          
  spline rho;            
  numero rV[Rho_imax];   
  numero phase;
  int i;
  Rho_density(rho, Z);                
  Rho_to_potential(Z,rho,rho,r,rV);   
  state.vint(r,rV,rho.n);
  state.dfree(E,phase,k);
  nmax=state.nmax;
  for(i=0; i<state.nmax; i++) {
    r[i]=state.rr[i];
    P[i]=state.P[i];  Q[i]=state.Q[i];
  }
  rho.free();  state.free();
}
#endif  
#endif  
#endif  
#endif  
#endif  
#ifndef jga_Mufpot        
#define jga_Mufpot 1      
                          
                          
class Mufpot {
public:
  int    init_flag;      
  int    read_flag;      
  numero d_error;        
  numero d_local;        
  Mufpot(void) {init_flag=read_flag=NR=0;  d_error=0.001;  d_local=5/a0_au;}
  
  void init(int n1, int n2, int *at, numero *x, numero *y, numero *z);
  void read(int n1, int n2, int *at, numero *x, numero *y, numero *z);
  void free();                                         
  
  int    N;              
  int    NR;             
  int    ngrid;          
  int    NX;             
                         
  numero *RX;            
                         
                         
  numero **Vcoul;        
  numero **rho;          
  numero **VH;           
  numero **VS;           
  numero **V;            
  numero Vmt0;           
                         
  int    *atZ;           
  int    *JRMT;          
  numero *RMT;           
  int    *NCON;          
  int    **IA;           
  int    **NA;           
  numero **AD;           
  
  void SUMAX(int ii, numero *pot_out, int IR);         
  numero pot_in(int i, int j, int k);                  
  void print(char *name, int cual);     
                                        
  
  void phases(int j,                                   
              numero Ei, numero Ef, int nE,            
              int lmax, numero *ph);                   
  void phases(int j, numero Ei, int lmax, numero *ph)
    {phases(j,Ei,Ei,1,lmax,ph);}
};
void Mufpot::free()
{
  int i;
  if(read_flag) {
    for(i=0; i<NR; i++) {
      if(init_flag)  delete [] V[i];
      delete [] IA[i];  delete [] NA[i];  delete [] AD[i];
    }
    delete [] NCON;  delete [] JRMT;
    delete [] IA;  delete [] NA;  delete [] AD;
    delete [] atZ;  delete [] RMT;  delete [] RX;
    if(init_flag)  delete [] V;
  }
  init_flag=read_flag=NR=0;
}
void Mufpot::print(char *name, int cual)
{
  int i,j;
  FILE *fout;
  if(!strcmp(name,"terminal") || !strcmp(name,"stdout"))  fout=stdout;
  else  fout=fopen(name, "w");
  for(i=0; i<78; i++)  fprintf(fout, "=");  fprintf(fout, "\n");
  if(read_flag && cual>=0) {
    fprintf(fout, "# atoms       = %3d\n", N);
    fprintf(fout, "# ineq. atoms = %3d\n", NR);
    for(j=0; j<NR; j++) {
      element_properties(atZ[j]);
      fprintf(fout, "***  %2d-type  ->  %2d shells   Rmt=%12.9f   %s\n",
                    j, NCON[j], RMT[j], element_symbol);
      for(i=0; i<NCON[j]; i++)
        fprintf(fout, "     shell=%2d   ta=%2d   n=%2d   d=%12.9f\n",
                    i, IA[j][i], NA[j][i], AD[j][i]);
  } }
  if(init_flag && cual<=0) {
    fprintf(fout, "Muffin-tin zero = %g eV\n", Vmt0*au_eV);
    fprintf(fout, " JRMT ->       ");
    for(i=0; i<NR; i++)  fprintf(fout, " %10d", JRMT[i]);
    fprintf(fout, "\n");
    for(j=0; j<NX; j++) {
      fprintf(fout, "%10.7f", RX[j]);
      for(i=0; i<NR; i++)  fprintf(fout, " %10.6f", RX[j]*V[i][j]);
      fprintf(fout, " // %3d\n", j);
  } }
  for(i=0; i<78; i++)  fprintf(fout, "=");  fprintf(fout, "\n");
  if(strcmp(name,"terminal") && !strcmp(name,"stdout"))  fclose(fout);
}
class atomic_environment {      
public:                         
  numero d;                     
  int    at, n;                 
  atomic_environment *next;     
  atomic_environment(void) {next=NULL;  n=0;}
  void free(void) {                                              
    if(next!=NULL) {next->free();  delete next;  next=NULL;}     
    n=0;
  }
  void add(numero &d_, int &at_, numero &d_error) {              
    if(at==at_ && ABS(d-d_)<=d*d_error)  n++;  else              
    if(next!=NULL)  next->add(d_,at_,d_error);  else {           
      next=new atomic_environment;                               
      next->d=d_;  next->at=at_;  next->n=1;  next->next=NULL;   
  } }
  int natoms(void) {                                             
    if(next!=NULL) return n+next->natoms();  return n;           
  }
  int find(atomic_environment &a, numero &d_error) {             
    if(n==a.n && at==a.at && ABS(a.d-d)<=d*d_error)  return 1;   
    if(next==NULL)  return 0;                                    
    return next->find(a, d_error);
  }
  void copy(atomic_environment &a) {                             
    free();                                                      
    n=a.n;  at=a.at;  d=a.d;
    if(a.next!=NULL) {
      next=new atomic_environment;
      next->copy(*(a.next));
    }  else  next=NULL;
  }
  void print(void) {                                             
    printf(" (%4.2f, %d, %d)", d, at, n);                        
    if(next!=NULL) next->print();                                
  }                                                              
};
int equal(atomic_environment &a, atomic_environment &b, numero &d_error)
{
  if(a.natoms()!=b.natoms())  return 0;                          
                                                                 
  atomic_environment *aa;  aa=&a;
  while(aa!=NULL) {if(b.find(*aa, d_error)==0) return 0;  aa=aa->next;}
  return 1;
}
void Mufpot::read(int n1, int n2, int *at, numero *x, numero *y, numero *z)
{
  int i,j,jj, eq, *at_, ZZ;   
  numero *RMT_;               
  numero xi,yi,zi, d,val;     
  atomic_environment *a;      
  atomic_environment aa,*pa;  
                              
  free();  read_flag=1;
  ngrid=250;                  
  RX=new numero [ngrid];
  for(i=0, d=-8.8; i<ngrid; i++) {RX[i]=exp(d);  d+=0.05;}
  a=new atomic_environment [n1];    for(i=0; i<n1; i++) a[i].next=NULL;
  at_=new int [n1];
  N=n1;  NR=0;
  for(i=0; i<n1; i++) {
    xi=x[i];  yi=y[i];  zi=z[i];                          
    aa.free();  aa.d=0;  aa.n=1;  aa.at=at[i];            
    for(j=0; j<n2; j++) if(j!=i) {                        
      d=sqrt(sqr(x[j]-xi)+sqr(y[j]-yi)+sqr(z[j]-zi));     
      if(d<d_local) aa.add(d, at[j], d_error);            
    }                                                     
    for(j=0, eq=0; j<NR && eq==0; j++)  {                 
      eq=equal(aa,a[j],d_error);
      if(eq==1) at_[i]=j;
    }
    if(eq==0) {a[NR].copy(aa);  at_[i]=NR;  NR++;}
  }
  for(i=0; i<n1; i++)  at[i]=at_[i];                      
  delete [] at_;                                          
  aa.free();
  NCON=new int [NR];  JRMT=new int [NR];
  IA=new int* [NR];  NA=new int* [NR];  AD=new numero* [NR];
  atZ=new int [NR];
  for(j=0; j<NR; j++) {
    pa=a+j;  NCON[j]=1;
    while(pa->next!=NULL) {NCON[j]++; pa=pa->next;}       
    IA[j]=new int [NCON[j]];
    NA[j]=new int [NCON[j]];
    AD[j]=new numero [NCON[j]];
    pa=a+j;  atZ[j]=pa->at;                               
    for(i=0; i<NCON[j]; i++) {                            
      IA[j][i]=pa->at;                                    
      NA[j][i]=pa->n;
      AD[j][i]=pa->d;
      pa=pa->next;
  } }
  for(j=0; j<NR; j++)  a[j].free();
  delete [] a;
  for(i=n1; i<n2; i++) {                                  
    jj=-1;  d=infinity;  zi=z[i];                         
    for(j=0; j<n1; j++)  if(atZ[at[j]]==at[i]) {          
      val=ABS(zi-z[j]);                                   
      if(val<d)  {d=val; jj=j;}                           
    }                                                     
    if(jj>=0)  at[i]=at[jj];  else  at[i]=-1;             
  }                                                       
  for(j=0; j<NR; j++)                                     
  for(i=1; i<NCON[j]; i++)                                
    if(IA[j][i]<1 || Rho_Zmax<IA[j][i])
      on_error("Mufpot::read", "atomic number out of range");
  RMT=new numero [NR];  RMT_=new numero [Rho_Zmax+1];     
  if(n1==1)  RMT[0]=Rho_radius[Rho_imax-1];  else  {
    for(j=1; j<=Rho_Zmax; j++)  RMT_[j]=-1;               
    for(j=0; j<NR; j++) {                                 
      RMT[j]=Rho_radius[Rho_imax-1];
      
      if(RMT_[atZ[j]]<0)  RMT_[atZ[j]]=Rho_rmt(atZ[j]);
      for(i=1; i<NCON[j]; i++) {
        ZZ=IA[j][i];
        if(RMT_[ZZ]<0)  RMT_[ZZ]=Rho_rmt(ZZ);             
        d=RMT_[atZ[j]]/RMT_[ZZ];                          
        d=AD[j][i] * d / (1+d);                           
        if(d<RMT[j])  RMT[j]=d;                           
  } } }
  delete [] RMT_;
}
void Mufpot::init(int n1, int n2, int *at, numero *x, numero *y, numero *z)
{
  spline rrho, ppot;
  numero Ze;
  int i,j,k, Zi;
  read(n1,n2, at, x,y,z);  init_flag=1;
  rho=new numero* [Rho_Zmax+1];  Vcoul=new numero* [Rho_Zmax+1];
  VH=new numero* [NR];  VS=new numero* [NR];  V=new numero* [NR];
  for(i=0; i<NR; i++) {                               
    VH[i]=new numero [ngrid];
    VS[i]=new numero [ngrid];
    V[i]=new numero [ngrid];
  }
  for(i=0; i<NR; i++)
  for(j=0; j<ngrid; j++) VH[i][j]=VS[i][j]=V[i][j]=0;
  for(i=1; i<=Rho_Zmax; i++) rho[i]=Vcoul[i]=NULL;
  for(i=0; i<NR; i++)                                 
  for(j=0; j<NCON[i]; j++)                            
    if(rho[IA[i][j]]==NULL) {
      Zi=IA[i][j];
      Rho_density(rrho, Zi);
      Ze=rrho.integ(0,rrho.b);
      rho[Zi]=new numero [ngrid];
      Vcoul[Zi]=new numero [ngrid];
      for(k=0; k<ngrid; k++) if(RX[k]<=rrho.b) {
        Vcoul[Zi][k]=(Poisson(RX[k],rrho)-Ze)/RX[k];  
        rho[Zi][k]=rrho.val(RX[k])/(4*pi*sqr(RX[k])); 
        NX=k;
      }  else  Vcoul[Zi][k]=rho[Zi][k]=0;
      rrho.free();
    }
  
  for(i=0; i<NR; i++) {
    JRMT[i]=int(20*(log(RMT[i])+8.8)+2); 
    SUMAX(0,VH[i],i);               
    SUMAX(1,VS[i],i);               
                                    
    for(j=0; j<NX; j++) VS[i][j]=Vxc_Barth_Hedin(VS[i][j],VS[i][j]/2);
  }
  
  
  
  
  for(i=0; i<NR; i++)
  for(j=0; j<JRMT[i]; j++)  V[i][j]=VH[i][j]+VS[i][j];
  NX=0;  for(i=0; i<NR; i++) if(NX<JRMT[i]) NX=JRMT[i];   
  for(i=1; i<=Rho_Zmax; i++) {
    if(rho[i]!=NULL)  delete [] rho[i];
    if(Vcoul[i]!=NULL)  delete [] Vcoul[i];
  }
  for(i=0; i<NR; i++) {delete [] VH[i];  delete [] VS[i];}
  delete [] VH;  delete [] VS;  delete [] rho;  delete [] Vcoul;
  Vmt0=0;  for(i=0; i<N; i++)  Vmt0+=V[at[i]][JRMT[at[i]]-1];  Vmt0=Vmt0/N;
  
  
  
  
  for(i=0; i<NR; i++)
  for(j=0; j<JRMT[i]; j++)  V[i][j]=V[i][j]-Vmt0;
}
numero Mufpot::pot_in(int i, int j, int ii)
{
  if(j>ngrid) return 0;
  if(ii==0)  return Vcoul[i][j-1];
  else       return rho[i][j-1];
}
void Mufpot::SUMAX(int ii, numero *pot_out, int IR)
{
  int ICON,JA,I,J,JTL,JBL;
  numero TOPX,XINT,ET,BLX,XBL,G,X,FZN,FZ1,FZ2,FZ3,TLX,XTL,C;
      ICON=IA[IR][0];
      TOPX=0.05*(NX-1)-8.8;
      for(J=1; J<=JRMT[IR]; J++)  pot_out[J-1]=pot_in(ICON,J,ii);
      if(NCON[IR]==1) goto gt16;
      JA=1;
gt2:  X=-8.8;  I=1;  ICON=IA[IR][JA];
gt3:  XINT=0;  ET=exp(X);  BLX=log(AD[IR][JA]-ET);
      if(BLX>=TOPX) goto gt12;
      JBL=int(2+20*(BLX+8.8));  XBL=0.05*(JBL-1)-8.8;  G=XBL-BLX;
      XINT+=0.5*G*pot_in(ICON,JBL,ii)*(2-20*G)*exp(2*XBL);
      XINT+=10*G*G*pot_in(ICON,JBL-1,ii)*exp(2*(XBL-0.05));
      TLX=log(AD[IR][JA]+ET);
      if(TLX<TOPX) goto gt6;
      JTL=NX;
      goto gt9;
gt6:  JTL=1+int(20*(TLX+8.8));
      if(JTL>=JBL) goto gt8;
      FZN=pot_in(ICON,JTL,ii)*exp(2*(XBL-0.05));
      FZ3=pot_in(ICON,JBL,ii)*exp(2*XBL);
      FZ2=FZN+20*(FZ3-FZN)*(TLX-XBL+0.05);
      FZ1=FZN+20*(FZ3-FZN)*(BLX-XBL+0.05);
      XINT=0.5*(FZ1+FZ2)*(TLX-BLX);
      goto gt12;
gt8:  XTL=0.05*(JTL-1)-8.8;  C=TLX-XTL;
      XINT+=0.5*C*pot_in(ICON,JTL,ii)*(2-20*C)*exp(2*XTL);
      XINT+=10*C*C*pot_in(ICON,JTL+1,ii)*exp(2*(XTL+0.05));
gt9:  if(JTL<=JBL) goto gt12;
gt10: XINT+=0.025*pot_in(ICON,JBL,ii)*exp(2*XBL);
      XINT+=0.025*pot_in(ICON,JBL+1,ii)*exp(2*(XBL+0.05));
      JBL++;
      if(JBL>=JTL) goto gt12;
      XBL+=0.05;
      goto gt10;
gt12: pot_out[I-1]+=0.5*XINT*NA[IR][JA]/(AD[IR][JA]*ET);
      if(I>=JRMT[IR]) goto gt14;
      I++;  X+=0.05;
      goto gt3;
gt14: if(JA+1>=NCON[IR]) goto gt16;
      JA++;
      goto gt2;
gt16: ;
}
void Mufpot::phases(int j,                           
              numero Ei, numero Ef, int nE,          
              int lmax, numero *ph)                  
{
  Salvat Dirac;
  c_au*=1000;
  Dirac.phases(RX,V[j],JRMT[j],atZ[j], Ei,Ef,nE, 0, lmax, ph);
  c_au/=1000;
}
#endif  
#ifdef jga_units          
#else                     
#define jga_units 1       
                          
                          
numero units_conv_fct(char *units)
{
  if(!strcmpC(units,"E(eV)"))    return 1/au_eV;     
  if(!strcmpC(units,"E(keV)"))   return 1000/au_eV;  
  if(!strcmpC(units,"E(au)"))    return 1.0;         
  if(!strcmpC(units,"E(Ry)"))    return 0.5;         
  if(!strcmpC(units,"k(A)"))     return a0_au;
  if(!strcmpC(units,"k(au)"))    return 1.0;
  if(!strcmpC(units,"l(A)"))     return 1/a0_au;
  if(!strcmpC(units,"l(au)"))    return 1.0;
  if(!strcmpC(units,"l(nm)"))    return 10/a0_au;
  if(!strcmpC(units,"l(mi)"))    return 10000/a0_au;
  if(!strcmpC(units,"v(au)"))    return 1.0;
  if(!strcmpC(units,"v(c)"))     return c_au;
  if(!strcmpC(units,"v(m/s)"))   return c_au/2.99792458e+8;
                                 return -1;
}
numero units_conversion_factor(char *units)
{
  numero val=units_conv_fct(units);
  if(val<0) on_error(foutput, "units_conversion_factor",
                              "wrong units specification", units);
  return val;
}
numero units_kvE(char *units, numero k, int i)
{
  k=k*units_conversion_factor(units);
  if(units[0]=='E' && i==0) {                   
    if(particle_type==electrones)  return sqrt(2*k);    
    if(particle_type==fotones)     return k/c_au;       
  }
  if(units[0]=='E' && i==1)                     
    return c_au*sqrt(1-1/sqr(k/sqr(c_au)+1));   
  if(units[0]=='l' && i==2)                     
    return 2*pi*c_au/k;                         
  if(units[0]!='E' && i==2)                     
    return sqr(k)/2;                            
  return k;
}
numero units_v(char *units, numero k) {return units_kvE(units,k,1);}  
numero units_E(char *units, numero k) {return units_kvE(units,k,2);}  
numero units_k(char *units, numero k) {return units_kvE(units,k,0);}  
char *str_units(char *units)
{
  if(!strcmpC(units,"E(eV)"))   return "eV";
  if(!strcmpC(units,"E(keV)"))  return "keV";
  if(!strcmpC(units,"E(au)"))   return "Hart.";
  if(!strcmpC(units,"E(Ry)"))   return "Rydb.";
  if(!strcmpC(units,"k(A)"))    return "1/A";
  if(!strcmpC(units,"k(au)"))   return "a.u.";
  if(!strcmpC(units,"l(A)"))    return "Angs.";
  if(!strcmpC(units,"l(nm)"))   return "nm";
  if(!strcmpC(units,"l(mi)"))   return "microns";
  if(!strcmpC(units,"l(au)"))   return "a.u.";
  if(!strcmpC(units,"v(au)"))   return "a.u.";
  if(!strcmpC(units,"v(c)"))    return "c";
  if(!strcmpC(units,"v(m/s)"))  return "m/s";
  on_error(foutput,"*str_units","wrong units specification",units);
  return units;
}
#endif  
#ifndef jga_algebra       
#define jga_algebra 1     
                          
                          
class algebra_variable {
 public:
  algebra_variable *next;
  char             *name;
  numero           a;
  algebra_variable(void) {name=NULL;}
} *algebra_var, *algebra_var0=NULL;
void algebra_free(void)
{
  while(algebra_var0!=NULL) {
    if(algebra_var0->name!=NULL)  delete [] algebra_var0->name;
    algebra_var=algebra_var0->next;  delete algebra_var0;
    algebra_var0=algebra_var;
} }
int algebra_find(char *name)
{
  algebra_var=algebra_var0;
  if(algebra_var==NULL)  return 0;
  while(strcmp(algebra_var->name,name)) {
    if(algebra_var->next==NULL)  return 0;
    algebra_var=algebra_var->next;
  }
  return 1;
}
#define algebra_error \
  on_error("algebraic evaluation", "wrong expression: ", ex);
int algebra_locate(char *ex, int i0, int i1, char c)
{
  int i=i1+1, p=0;  
  int cond;         
  do {
    i--;
    if(ex[i]=='(' || ex[i]=='[' || ex[i]=='{') p--;
    if(ex[i]==')' || ex[i]==']' || ex[i]=='}') p++;
    if(p<0) algebra_error
    cond=(ex[i]==c);
    if(cond) if(c=='+' || c=='-') {
      if(i>i0) if(ex[i-1]=='e' || ex[i-1]=='d')
      if(i>i0+1) if(ex[i-2]=='.' || ('0'<=ex[i-2] && ex[i-2]<='9')) cond=0;
    }
  } while((!cond || p!=0)  && i>i0);
  if(p!=0) algebra_error
  if(!cond) return -1;
  return i;
}
int algebra_locate(char *ex, char c)
  {return algebra_locate(ex,0,strlen(ex)-1,c);}
int algebra_locate_absolute(char *ex, int i0, int i1, char c)
{
  int i=i0-1;  do i++; while(ex[i]!=c && i<i1);
  if(ex[i]!=c) return -1;
  return i;
}
numero algebra_eval(char *ex, int i0, int i1)
{
  char c,d;  int i;  if(i0>i1) return 0;  c=ex[i0];  d=ex[i1];
  numero val;
  if(c=='*' || c==':' || c=='/' || c=='^' || c==')' || c==']' || c=='}' ||
     d=='*' || d==':' || d=='/' || d=='^' || d=='(' || d=='[' || d=='{')
    algebra_error
  if((i=algebra_locate(ex,i0,i1,'+'))>=0)
    return algebra_eval(ex,i0,i-1)+algebra_eval(ex,i+1,i1);
  if((i=algebra_locate(ex,i0,i1,'-'))>=0)
    return algebra_eval(ex,i0,i-1)-algebra_eval(ex,i+1,i1);
  if((i=algebra_locate(ex,i0,i1,'*'))>=0)
    return algebra_eval(ex,i0,i-1)*algebra_eval(ex,i+1,i1);
  if((i=algebra_locate(ex,i0,i1,':'))>=0)
    return algebra_eval(ex,i0,i-1)*algebra_eval(ex,i+1,i1);
  if((i=algebra_locate(ex,i0,i1,'/'))>=0)
    return algebra_eval(ex,i0,i-1)/algebra_eval(ex,i+1,i1);
  if((i=algebra_locate(ex,i0,i1,'^'))>=0)
    return pow(algebra_eval(ex,i0,i-1),algebra_eval(ex,i+1,i1));
  if(c=='(' && d==')')  return algebra_eval(ex,i0+1,i1-1);
  if(c=='[' && d==']')  return algebra_eval(ex,i0+1,i1-1);
  if(c=='{' && d=='}')  return algebra_eval(ex,i0+1,i1-1);
  char *aux;  aux=new char [i1-i0+2];
  for(i=i0; i<=i1; i++) aux[i-i0]=ex[i];  aux[i1-i0+1]=0;  i1=i1-i0;  i0=0;
  if(('a'<=c && c<='z') || ('A'<=c && c<='Z')) {d=aux[i1];
    if(d==')' || d==']' || d=='}') {
      if((i=algebra_locate_absolute(aux,i0,i1,'('))<0)
      if((i=algebra_locate_absolute(aux,i0,i1,'['))<0)
      if((i=algebra_locate_absolute(aux,i0,i1,'{'))<0)
        algebra_error
      val=algebra_eval(aux,i+1,i1-1);
      aux[i]=0;
      #ifdef jga_bessel
      if(!strcmp(aux,"J0"))    val=besselJJ(0,val);  else
      if(!strcmp(aux,"J1"))    val=besselJJ(1,val);  else
      if(!strcmp(aux,"J2"))    val=besselJJ(2,val);  else
      if(!strcmp(aux,"J3"))    val=besselJJ(3,val);  else
      if(!strcmp(aux,"Y0"))    val=besselYY(0,val);  else
      if(!strcmp(aux,"Y1"))    val=besselYY(1,val);  else
      if(!strcmp(aux,"Y2"))    val=besselYY(2,val);  else
      if(!strcmp(aux,"Y3"))    val=besselYY(3,val);  else
      if(!strcmp(aux,"j0"))    val=besselj(0,val);   else
      if(!strcmp(aux,"j1"))    val=besselj(1,val);   else
      if(!strcmp(aux,"j2"))    val=besselj(2,val);   else
      if(!strcmp(aux,"j3"))    val=besselj(3,val);   else
      if(!strcmp(aux,"y0"))    val=bessely(0,val);   else
      if(!strcmp(aux,"y1"))    val=bessely(1,val);   else
      if(!strcmp(aux,"y2"))    val=bessely(2,val);   else
      if(!strcmp(aux,"y3"))    val=bessely(3,val);   else
      #endif
      #ifdef jga_legendre
      if(!strcmp(aux,"P00"))   val=legendre(0,0,val);  else
      if(!strcmp(aux,"P10"))   val=legendre(1,0,val);  else
      if(!strcmp(aux,"P11"))   val=legendre(1,1,val);  else
      if(!strcmp(aux,"P20"))   val=legendre(2,0,val);  else
      if(!strcmp(aux,"P21"))   val=legendre(2,1,val);  else
      if(!strcmp(aux,"P22"))   val=legendre(2,2,val);  else
      #endif
      if(!strcmp(aux,"ln"))    val=log(val);    else
      if(!strcmp(aux,"log"))   val=log(val);    else
      if(!strcmp(aux,"exp"))   val=exp(val);    else
      if(!strcmp(aux,"sin"))   val=sin(val);    else
      if(!strcmp(aux,"cos"))   val=cos(val);    else
      if(!strcmp(aux,"tan"))   val=tan(val);    else
      if(!strcmp(aux,"asin"))  val=asin(val);   else
      if(!strcmp(aux,"acos"))  val=acos(val);   else
      if(!strcmp(aux,"atan"))  val=atan(val);   else
      if(!strcmp(aux,"abs"))   val=ABS(val);    else
      if(!strcmp(aux,"mod"))   val=ABS(val);    else
      if(!strcmp(aux,"ABS"))   val=ABS(val);    else
      if(!strcmp(aux,"sqr"))   val=sqr(val);    else
      if(!strcmp(aux,"sqrt"))  val=sqrt(val);   else
      if(!strcmp(aux,"sinh"))  val=(exp(val)-exp(-val))/2;   else
      if(!strcmp(aux,"cosh"))  val=(exp(val)+exp(-val))/2;   else
      on_error("algebraic evaluation", "undefined function", aux);
    } else
    if(!strcmp(aux,"pi")) val=pi; else
    if(!strcmp(aux,"PI")) val=pi; else
    if(!strcmp(aux,"Pi")) val=pi; else
    if(!strcmp(aux,"e"))  val=exp(1.0); else
    if(!strcmp(aux,"E"))  val=exp(1.0); else
    if(!strcmp(aux,"eV"))  val=1/au_eV; else
    if(!strcmp(aux,"nm"))  val=nm; else
    if(!strcmp(aux,"c"))   val=c_au; else
    if(!strcmp(aux,"infinite")) val=infinity; else
    if(!strcmp(aux,"infinity")) val=infinity; else
    if(!strcmp(aux,"zero"))     val=0;        else
    if(algebra_find(aux)==0) on_error("algebraic evaluation",
       "undefined variable", aux);
    else  val=algebra_var->a;
    delete [] aux;  return val;
  }
  val=read_numero(aux);
  delete [] aux;  return val;
}
numero algebra_eval(char *ex)
{
  int i,l=strlen(ex);  char c;
  for(i=0; i<l; i++) {
    c=ex[0];
    if(c!='+' && c!='-' && c!='*' && c!=':' && c!='/' && c!='^' && c!='.' &&
       c!='(' && c!=')' && c!='[' && c!=']' && c!='{' && c!='}' &&
       (c<'a' || c>'z') && (c<'A' || c>'Z') && (c<'0' || c>'9'))
      algebra_error
  }
  return algebra_eval(ex,0,l-1);
}
int algebra_verbose=0;
void algebra_define(char *name, char *ex)
{
  numero val=algebra_eval(ex);
  if(algebra_find(name)==0) {
    algebra_var=new algebra_variable;
    algebra_var->name=new char [strlen(name)+1];
    strcpy(algebra_var->name,name);
    algebra_var->next=algebra_var0;
    algebra_var0=algebra_var;
  }
  algebra_var->a=val;
  if(algebra_verbose)
    printf(">>> %s = %s = %g <<<\n", algebra_var->name,ex,algebra_var->a);
}
numero alread_numero(char *name)
{
  strcpy(read_command,name);
  return algebra_eval(read_command);
}
int alread_int(char *name)
{
  strcpy(read_command,name);
  return ((int) algebra_eval(read_command));
}
numero alread_numero(FILE *f)
{
  read_name(f,read_command);
  return algebra_eval(read_command);
}
int alread_int(FILE *f)
{
  read_name(f,read_command);
  return ((int) algebra_eval(read_command));
}
#endif  
class green;               
class green_z;             
class rotation_th;         
class rotation_fi;         
class t_matrix;            
class vector;              
class cluster;             
class atomic_scattering;   
class propagation;         
class multiple_scattering; 
class final_state;         
class calculation;         
class electron_bin;        
class green {        
public:
  numero *re,*im;    
  numero d,th,fi;    
  green(void) {d=th=fi=0;}
  green *next;       
  void free(void) {delete [] re;  delete [] im;}
  void init(complex k, numero dd, numero tth, numero ffi, int periodic);
  void prod(vector &a);
};
class green_z {      
public:
  numero *re,*im,d;  
  int trans;         
  green_z(void);
  green_z *next;     
  void free(void) {delete [] re;  delete [] im;}
  void init(complex k, numero dd, int translation_type, int att_type);
};
class rotation_th {  
public:
  numero *val,th;    
  rotation_th *next; 
  rotation_th(void);                 
  void free(void) {delete [] val;}
  int init(numero tth);  
};
class rotation_fi {      
public:
  numero *re,*im,fi;     
  rotation_fi *next;     
  rotation_fi(void);
  void free(void) {delete [] re;  delete [] im;}
  int init(numero ffi);  
};
class t_matrix {         
public:
  numero *re,*im;        
  numero rmt;            
  int lmax;              
  int n_c;               
  int ta;                
  int select;            
  t_matrix *next;        
  void  free(void) {if(n_c>0) {delete [] re;  delete [] im;  n_c=0;}}
  t_matrix (void) {n_c=0;  ta=-1;  lmax=-1;  rmt=0;}
};
class vector {       
public:
  int lmax;          
  int n_c;           
  numero *re,*im;    
  vector(void);
  void free(void);
  void alloc(int lmax_, int n_c_);      
  void copy(vector &a);
  void copy_to_array(complex *a);
  void copy_to_array(complex *a, int cc);
  void copy_to_array(complex *a, int lmax_, int cc);
  void copy_from_array(complex *a, int lmax_, int n_c_);
  void substitute_lmax(int lmax_);
  void substitute_n_c(int n_c_, int c_flag);
  int spin_rotation(numero al, numero be, numero ga);
  friend vector operator*(numero x, vector &b);
  int add(vector &b);
  int ext(numero x);
  int Rth(rotation_th *R, int op);
  int Rfi(rotation_fi *R, int op);
  int Gz(green_z *G, int lmax_final);
  int tmat(t_matrix *t);
private:
  void copy(vector &a, int lmax_, int n_c_, int c_flag);
  void tmat_diagonal(t_matrix *t);
  void tmat_full(t_matrix *t);
  void tmat_so(t_matrix *t);
};
class cluster {
public:
  
  int    n,                    
         nn,                   
         nnn;                  
  numero *x,*y,*z;             
  int    *atom;                
  int    *atZ;                 
                               
  
  numero Rmax;                 
  numero x0,y0,z0;             
  
  int    oriented;             
  numero *th,*fi;              
  
  numero eps_d;                
  numero eps_th;               
  
  int    n_emit,               
         *emit,                
         *emit_s;              
  void add_emitter(int a);                                   
  void add_emitter(int a, int sc);
  void add_emitter_plane(numero x_, numero y_, numero z_,    
                         numero th_, numero fi_);            
  void add_emitter_column(numero x_, numero y_);             
  int is_emitter(int i);       
  void free_emitters(void);    
  void all_emitters(int ta, int sc);   
  void report_emitters(void);  
  
  numero zs;                   
  int    surface;              
  
  cluster(void);               
  void free(void);             
  void free_coordinates();     
  int set_n(void);             
  int  read_atomic_type(FILE *fin);                   
  void read(FILE *fin, char *my_file);                
  void write(char *name, char *units);                
  void read(FILE *fcoordinates);
  void write(FILE *fout, char *units);
  numero distance(numero x0_, numero y0_, numero z0_, 
                  numero x1, numero y1, numero z1);   
  int get_atom(numero x_, numero y_, numero z_);      
  int get_atom(numero x_, numero y_, numero z_, numero err);
  void sort(int emitter);                             
  int add_atom(numero x_, numero y_, numero z_,       
               numero th_, numero fi_, int at);
  int delete_atom(numero x_, numero y_, numero z_);   
  void add_layer(numero x_, numero y_, numero z_,     
                 numero a, numero b,
                 numero al, numero be,
                 numero th_, numero fi_, int at);
  void delete_layer(numero x_, numero y_, numero z_,  
                    numero a, numero b,
                    numero al, numero be);
  void add_surface(numero x_, numero y_, numero z_,   
                 numero a, numero c, char *name,
                 numero th_, numero fi_, int at);
  void add_half_layer(numero x_, numero y_, numero z_,
                      numero a, numero b,
                      numero al, numero be, numero ff,
                      numero th_, numero fi_, int at);
  void add_row(numero x_, numero y_, numero z_, numero a,
                numero tt, numero ff, numero th_, numero fi_, int at);
  void delete_row(numero x_, numero y_, numero z_,
                   numero a, numero tt, numero ff);
  void displace(numero xx, numero yy, numero zz);     
  void rotate(int axis, numero ang);                  
  void rotate(double &xx, double &yy, double ang);    
  int delete_non_emitters(void);                      
};
class atomic_scattering {
public:
  atomic_scattering(void);
  
  atomic_scattering *next;   
                             
                             
  int ta;                    
                             
                             
                             
  int     select;            
  int     n_c;               
  int     lmax;              
                             
                             
  numero  *re, *im;          
  int magnetic;              
                             
                             
  int     n_k;               
  numero  *k;                
                             
  numero  mass;              
  
  void free(void);                     
  void init_t_matrix(t_matrix *t);     
  int read(FILE *fin, char *name);     
  int scan(FILE *fout, int sel);       
  
  int comp;                              
  spline *rr, *ii;                       
  void read_phases(FILE *fin, char *units);
  complex *tl_phases(int n_c_, int &lmax_);
  void cut_tl_phases(complex *tl, int n_c_, int &lmax_);
  numero Gaunt_j(int l, int j, int mu, int lp, int jp, int lpp);
  complex *tl_temperature(complex *tl, int n_c_, int &lmax_);
  complex *tl_temperature_so(complex *tl, int &lmax_);
  numero tl_deltal(complex tl) {return atan2(imag(tl), real(tl));}
  void init_t_matrix_diag(t_matrix *t);
  void init_t_matrix_diag(complex *tl, t_matrix *t);
  void init_t_matrix_so(t_matrix *t);
  void init_t_matrix_so(complex *tl, t_matrix *t);
  
  numero a,V;                            
  complex *tl_sphere(int lmax_);
  void init_t_matrix_sphere(t_matrix *t);
  
  int    n_rad;                          
  numero *rad;                           
  int    atZ;                            
  complex *tl_potential(int n_c_, int &lmax_);   
  void init_t_matrix_potential(t_matrix *t);
  void read_potential(FILE *fin);
  
  void read_full(FILE *fin, char *units);
  void init_t_matrix_full(numero kk, t_matrix *t);
  void free_full(void);
  
  numero *copy_(numero *a_, numero *b_, int n, numero x);
};
#define n_Debye_temperature_elements 100
class propagation {
public:
  
  int     n_d,                
          n_th,               
          n_fi,               
          n_t;                
  numero  dmax;               
  complex k;                  
  
  int     translation_type;   
  int     attenuation_type;   
  numero  T;                  
  numero  T_Debye[n_Debye_temperature_elements];
                              
                              
  numero  V0;                 
  numero  refraction;         
  spline  iimfp;              
  int     iimfp_flag;         
  numero  Vi;                 
  void    scan_imfp(char *name);      
  void    scan_imfp(FILE *fout);      
  numero  iimfp_TPP(numero kr);
  numero  TPP_rho, TPP_Nv, TPP_Ep, TPP_Eg;
  numero  screening_length;   
  int     scattering_so;      
  void set_k(int ik);         
  void set_k(complex k_);
  void read_imfp(FILE *fin, char *my_file);
  void set_muffin_tin_potentials(void);
  
  green_z     *Gz,*Gz0;       
  rotation_th *Rth,*Rth0;     
  rotation_fi *Rfi,*Rfi0;     
  t_matrix    *t_mat,
              *t_mat0;        
  numero      d,              
              th,             
              fi;             
  int         ifi;            
  
  int n_a;                       
  atomic_scattering *a, *a0;     
  numero epsilon;                
  
  propagation(void);
  void R_1TR(vector &gg, int i, t_matrix *tt);
  void R_1TR(vector &gg, int i);
  void Rab(vector &gg);
  void Rab_1(vector &gg);
  void propagate(vector &gg, int lmax_final);
  void free_Green_functions(void);
  void free_rotations(void);
  void free_t_matrix(void);
  void report_links(void);
  int get_ij(int i, int j);
  void get_ij(void);
  int get_d(void);
  int get_th(numero th_);
  int get_fi(numero fi_);
  int get_t_matrix(int ta);
  int get_scatterer(int ta, int opt);
  int get_scatterer(int ta) {return get_scatterer(ta,1);}
  int free_scatterer(int ta);
  void free_scatterers(void);
  void add_scatterer(int ta);
  void add_screened_scatterer(int ta, int ta_old);
  void scan_scatterer(int ta, char *name, int sel);
  void set_lmax(void);        
  int get_n_c(void);          
  
  int dim;                    
  matrix aa;                  
  int lmax_aa;                
  void get_exact(int n_c, int lmax_);  
  void free_exact();
  char aa_output[100];
  int  aa_output_flag;
  
  int    separation;          
  int    n_G;                 
  green  *Gab,*Gab0;
  numero *Gab_cg;
  int    Gab_cg_flag;
  int get_Gab(void);
  int    nn_d,nn_th,nn_fi,nn_t,nn_G,nn_c;  
};
                  
class wave_function {         
public:                       
                              
  int n;                      
  vector *g;                  
  wave_function(void) {n=0;}  
  void alloc(int n_) {                                       
     free();  n=n_;  g=new vector [n];
     int i;  for(i=0; i<n; i++) {g[i].lmax=-1;  g[i].n_c=1;}
  }
  void clear(void) {int a; for(a=0; a<n; a++) g[a].free();}  
  void free(void) {clear();  if(n>0) delete [] g;  n=0;}     
};
                  
class multiple_scattering {
public:
  
  wave_function o,n;  
  
  int     method;     
  numero  eta;        
  
  multiple_scattering(void);
  void next_scattering(void);
  void exact(void);
  void free(void);
  void clear(void);
  void alloc(void);
  
  void write(FILE *fout, int o_, int max_order_);
  int read(FILE *fin, int o_);
  int max_order;      
};
                  
class array {            
public:
  complex *a;            
  int    n;              
  array(void) {n=0;}
  void free(void) {if(n>0) delete [] a;  n=0;}
  void alloc(int n_) {free();  n=n_;  a=new complex [n];
                      int i;  for(i=0; i<n; i++) a[i]=0;}
};
class final_state {
public:
  
  int     n_th;             
  int     n_fi;             
  numero  *th, *fi;         
  
  numero  *th_out,          
          *transmission;    
  
  final_state(void);
  void free(void);
  void init_th(numero thi, numero thf, int nth);
  void init_phi(numero fii, numero fif, int nfi);
  void init_refraction(
            numero refraction);  
  void init_transmission(        
            numero refraction);
  
  void free_mesh(void);
  void init_mesh(int n_1_, int n_2_);
  void init_mesh(numero x0, numero y0, numero z0,
                 numero x1, numero y1, numero z1,
                 numero x2, numero y2, numero z2, int n_1_, int n_2_);
  void init_mesh(numero r,
                 numero th1, numero th2, numero fi1, numero fi2,
                 int n_1_, int n_2_);
  int mesh_flag;
  numero *x,*y,*z;          
  int    n_1,n_2;           
  
  void init_Ylm(void);
  void init_Ylm(numero th_, numero fi_, complex *Ylm_);
  complex *Ylm(int i, int j); 
  vector Ylm(complex *Ylm_);  
  void free_Ylm(void);
  array  *Ylm0_th;          
  array  *Ylm0_fi;          
  int    Ylm0_th_flag;      
  int    Ylm0_fi_flag;      
};
                  
class calculation {
public:
  
  int         n_final,   
              n_ang,     
              n_mesh,    
              n_initial, 
              n_initial_used,  
              n_c;       
  numero  thave;         
  numero  *III0;
  numero  IIIthfi(int no, numero th, numero fi);
  numero  IIIave(int no, numero th, numero fi);
  int         n_atoms_last;    
  numero  e_th, e_fi,    
          e_tth, e_ffi,  
          alpha, delta;  
  complex eps[3];        
  int     integ;         
                         
                         
  numero  *k, *kk;       
  char    *k_units;      
  numero  V0;            
  int     n_k;           
  int     k_flag;        
  complex *kc;           
  numero  beta;          
  int     mode_beta,     
          mode_sample,
          mode_eps;
  int     beamline;      
  numero  beamline_th, beamline_fi;  
  
  int     *orders,       
          n_orders,      
          max_order;     
  int n_orders_eff(void);    
  
  char    *wf_in,        
          *wf_out;       
  FILE    *wf_fin,       
          *wf_fout;      
  int     ms;            
  
  array   *results;      
  FILE    *fout;         
  int     verbose;       
  int     select_wf;     
  int     n_intensities; 
  int     n_f;           
  
  calculation(void);
  void free_results(void);
  void free_k(void);
  void free_orders(void);
  void free(void);
  void alloc(int n_intial_);
  
  void init_k(numero kki, numero kkf, int n_k_, char *k_units_);
  void init_k(complex *kc_, int n_k_, char *k_units_);
  void init_incidence(numero e_th_, numero e_fi_);
  void init_polarization(numero alpha_, numero delta_);
  void init_orders(int n_orders_, int *o);
  void init_final_dist(void);
  void init_final_f(int n_f_);
  void init_mode_sample(int mode_sample_);
  void init_mode_beta(int mode_beta_);
  void init_beamline(int beamline_);
  void init(char *name, char *title);
  void init(char *name);
  
  int lmax;                     
  int lmax_inter;               
  int dim;                      
  int center;                   
  void init_wf_tmat(int c, int l, int m);
  complex *wave_function_tmat(complex *wf);
  void scan_tmat(char *name, int lmax_);
  
  complex *wave_function(complex *wf);
  int evaluation(int initial);
  void eval_recursion(int initial);
  
  complex *wave_function(int i, int j, complex *wf);
  void intensity(numero *intensity_, complex *wf);
  
  void write(int ik);
  
  void set_polarization(void);
  void set_polarization(int pol);
  void eval_polarization(void);
private:
  int  polarization_flag;   
  void write_title(FILE *fout, char *title);
  void write_orders(FILE *fout);
  void write_ang(FILE *fout, char *title);          
  void write_mesh(FILE *fout, char *title);
  void write_f(FILE *fout, char *title);
  void write_k(FILE *fout_, int ik);                
  void write_ang(FILE *fout, int ik);               
  void write_mesh(FILE *fout, int ik);
  void write_f(FILE *fout, int ik);
  void write_numeros(FILE *fout_, numero *intensity_, complex *wf);
  void write_numeros(FILE *fout_, numero val);
  void intensity(numero *intensity_, int no, int i, int j, complex *wf);
};
class electron_bin {
public:
  
  int     spin_flag;          
  numero  th_spin,            
          fi_spin;            
  
  int     molecule;           
  
  int     n0,l0,              
          m0,                 
          ms0,                
          j0, mj0,            
          ljmj;               
  
  int     n_c;                
  int     excitation;         
  int     n_rmat;             
  splinec *rmat;              
  int     emitter;            
  
  int    CR_nterms, CR_flag;
  numero CR_C[10], CR_xhi[10];
  int    CR_n[10];
  
  numero intensity(complex *wf);
  void wf_far(int i, int j, complex *wf);
  void wf_far(numero th, numero fi, complex *Ylm, complex *wf);
  void wf_far(vector &a, complex *Ylm, complex *wf);
  void wf_near(int i, int j, complex *wf);
  
  int scan_PD(char *name);           
  int scan_MS(char *name, int nwf);  
  int wf_ns;                         
  
  complex  rrm[4];                   
  void init_wf_PD(int jm, int mjms, complex *eps);   
  void init_wf_lmms(int m0_, int ms0_, complex *eps);
  void init_wf_ljmj(int j0_, int mj0_, complex *eps);
  complex ang_mat_el(int m0_, int l, int m, complex *eps);   
  void read_rmat(FILE *fin);      
  void free_rmat(void);
  complex rmat_el_MT(int l);      
  complex rmat_el_CR(int l);      
  numero  CR_wf(numero r);        
  spline  P0,Q0;                  
  void init_rmat_MT(void);
  void free_rmat_MT(void)  {P0.free();  Q0.free();}
  complex *manual;                
  int manual_flag;                
  int lmax_manual;                
  int n_c_manual;                 
  void init_wf_manual(void);      
  
  int      leed_i, leed_j;    
  numero   leed_th;           
  numero   leed_x, leed_y;    
  numero   leed_a, leed_b;    
  numero   ax,ay,bx,by;       
  numero   Gax,Gay,Gbx,Gby;   
  numero   area;              
  numero   k0x,k0y,k0z;       
  void surface_lattice(numero a, numero b, numero al, numero be);
  void init_leed_emitters(void);
  int init_wf_LEED(int spin_);        
  int scan_DLEED(char *name);         
  int scan_LEED(char *name);          
  void wf_LEED(complex *wf);          
  
  void program(char *my_file);
  
  electron_bin(void)
  {
    th_spin=fi_spin=0;  spin_flag=0;
    ljmj=0;  n0=1; l0=0;
    m0=infinite_int+1;       
    ms0=1;
    j0=infinite_int+1;  mj0=infinite_int+1;    
    n_c=1;  wf_ns=1;
    manual=NULL;  manual_flag=0;
    rmat=NULL;
    leed_i=leed_j=0;
    leed_x=leed_y=0;
    surface_lattice(5/a0_au,5/a0_au, 0,pi/2);
    CR_flag=0;               
    molecule=0;
} };
void init_lmax(int &lmax);       
void init_lmax(void);            
void free_lmax(void);
void report_size(void);
cluster              coord;         
propagation          scat;          
multiple_scattering  iter;          
final_state          final;         
calculation          calc;          
electron_bin         electron;      
int size_of_rotation_th=0,  
    size_of_rotation_fi=0,  
    size_of_green_z=0,      
    size_of_lmax=-1;        
int interpolation_type=0;  
#define str_length_files     100
#define str_length_comments  100
#define str_length_units      10
int lmax_eff=-1,            
    lmax_max=15;            
numero *Gaunt_integrals,    
       *cgcg12;             
int    *lp12;               
void init_lmax(int &lmax)
{
  free_lmax();   size_of_lmax=lmax_eff=lmax=(lmax>lmax_max)?lmax_max:lmax;
  init_lmax();
}
void init_lmax(void)
{
  int l1,l2,l3,m1,m2,p, n_CG;
  for(l1=0, size_of_green_z=0, n_CG=0; l1<=lmax_eff; l1++)
  for(m1=-l1; m1<=l1; m1++)
  for(l2=(m1>0?m1:-m1); l2<=lmax_eff; l2++, size_of_green_z++)
  for(l3=(l1>l2?l1-l2:-l1+l2); l3<=l1+l2; l3+=2) n_CG++;
  Gaunt_integrals = new numero [n_CG];
  numero pi_factor=sqrt(4*pi);  int i;
  for(l1=0, i=0; l1<=lmax_eff; l1++)                 
  for(m1=-l1; m1<=l1; m1++)                          
  for(l2=(m1>0?m1:-m1); l2<=lmax_eff; l2++)          
  for(l3=(l1>l2?l1-l2:-l1+l2); l3<=l1+l2; l3+=2)
    Gaunt_integrals[i++]=pi_factor*sqrt(2*l3+1)*Gaunt(l1,l3,l2,m1,0,m1);
  for(l1=0, size_of_rotation_th=0; l1<=lmax_eff; l1++)
  for(m1=-l1; m1<=l1; m1++)
  for(m2=-l1; m2<=l1; m2++) size_of_rotation_th++;   
                                                     
                                                     
  size_of_rotation_fi=2*lmax_eff+1;
  lp12 = new int [lmax_eff+3];
  for(l1=0; l1<=lmax_eff+2; l1++)  lp12[l1]=(l1+1)*(l1+1);
  cgcg12 = new numero [8*lp12[lmax_eff]];
  for(l1=p=0; l1<=lmax_eff; l1++)
  for(m1=-l1; m1<=l1; m1++, p+=8)  {
    l2=l1+l1;  m2=m1+m1;
    cgcg12[p  ]=CGspin1_2(l2,l2+1, 1,m2+1) * CGspin1_2(l2,l2+1, 1,m2+1); 
    cgcg12[p+1]=CGspin1_2(l2,l2-1, 1,m2+1) * CGspin1_2(l2,l2-1, 1,m2+1); 
    cgcg12[p+2]=CGspin1_2(l2,l2+1,-1,m2-1) * CGspin1_2(l2,l2+1,-1,m2-1); 
    cgcg12[p+3]=CGspin1_2(l2,l2-1,-1,m2-1) * CGspin1_2(l2,l2-1,-1,m2-1); 
    cgcg12[p+4]=CGspin1_2(l2,l2+1, 1,m2+1) * CGspin1_2(l2,l2+1,-1,m2+1); 
    cgcg12[p+5]=CGspin1_2(l2,l2-1, 1,m2+1) * CGspin1_2(l2,l2-1,-1,m2+1); 
    cgcg12[p+6]=CGspin1_2(l2,l2+1,-1,m2-1) * CGspin1_2(l2,l2+1, 1,m2-1); 
    cgcg12[p+7]=CGspin1_2(l2,l2-1,-1,m2-1) * CGspin1_2(l2,l2-1, 1,m2-1); 
} }
void free_lmax(void)
{
  if(lmax_eff>=0) {
    delete [] Gaunt_integrals;  delete [] lp12;
    delete [] cgcg12;
    lmax_eff=-1;
} }
void report_size(void)
{
  fprintf(foutput,"%s%s%4d%s%4d%s%4d%s%4d\n",
            "--- size of different elements",
          "\n              lmax              = ", size_of_lmax,
          "\n              Green's function  = ", size_of_green_z,
          "\n              Rotation (0,th,0) = ", size_of_rotation_th,
          "\n              Rotation (0,0,fi) = ", size_of_rotation_fi);
}
calculation::calculation(void)
{                                       
  n_final=n_orders=n_c=n_k=0;
  thave=0;                              
  init_incidence(0,0);                  
  init_polarization(0,0);               
  beta=0;                               
  mode_beta=0;                          
  mode_sample=0;                        
  mode_eps=1;
  fout=NULL;                            
  verbose=1;                            
  k_units=new char [str_length_units];
  integ=0;                              
  n_atoms_last=0;                       
  select_wf=0;                          
  n_intensities=1;                      
  n_ang=n_mesh=n_f=0;
  n_initial=n_initial_used=0;
  wf_in=new char [str_length_files];    strcpy(wf_in,"none");
  wf_out=new char [str_length_files];   strcpy(wf_out,"none");
  wf_fin=NULL;  wf_fout=NULL;
  ms=0;                                 
  beamline=0;                           
}
void calculation::init(char *name)
{
  close_file(fout);
  close_file(wf_fin);
  close_file(wf_fout);
  wf_fin =open_file(foutput,wf_in, "r");    
  wf_fout=open_file(foutput,wf_out,"a");    
  fout   =open_file(foutput,name,  "w");    
}
void calculation::init(char *name, char *title)
{
  init(name);
  if(n_ang>0) {
    if(verbose)    write_ang(foutput,title);        
    if(fout!=NULL) write_ang(fout,title);
  } else
  if(n_mesh>0) {
    if(verbose)    write_mesh(foutput,title);
    if(fout!=NULL) write_mesh(fout,title);
  } else
  if(n_f>0) {
    if(verbose)    write_f(foutput,title);
    if(fout!=NULL) write_f(fout,title);
} }
void calculation::init_final_dist(void)
{
  n_f=0;
  n_ang=final.n_th * final.n_fi;
  if(n_ang>0)  n_mesh=0;
  else         n_mesh=final.n_1 * final.n_2;
  n_final=n_ang+n_mesh;
}
void calculation::init_final_f(int n_f_)
{
  n_ang=n_mesh=0;  n_final=n_f=n_f_;
}
void calculation::init_orders(int n_orders_, int *o)
{
  int i,j,l,m;
  free_orders();
  n_orders=n_orders_;
  orders=new int [n_orders];               
  for(i=n_orders-1; i>=0; i--) {           
    l=-1;
    for(j=0; j<n_orders; j++)  if(l<o[j]) {l=o[j]; m=j;}
    orders[i]=l;  o[m]=-1;
  }
  delete [] o;
  max_order=orders[n_orders-1];
}
void calculation::init_k(numero kki, numero kkf, int n_k_, char *k_units_)
{
  numero dkk;
  int i;
  free_k();
  if(n_k_<1) n_k=1;  n_k=n_k_;  k_flag=1;
  strcpy(k_units,k_units_);
  kk=new numero [n_k];  k=new numero [n_k];
  if(n_k>1) dkk=(kkf-kki)/(n_k-1);  else dkk=0;
  for(i=0; i<n_k; i++) {
    kk[i]=kki+i*dkk;
    k[i]=units_k(k_units,kk[i]);
    if(k[i]<=0)  on_warning(foutput, "calculation::init_k", "momentum k<0");
} }
void calculation::init_k(complex *kc_, int n_k_, char *k_units_)
{
  free_k();
  n_k=n_k_;
  strcpy(k_units,k_units_);
  kc=kc_;  k_flag=2;
}
void calculation::init_incidence(numero e_th_, numero e_fi_)
{
  normalize_angles(e_th_,e_fi_);   e_th=e_tth=e_th_;   e_fi=e_ffi=e_fi_;
}
void calculation::init_polarization(numero alpha_, numero delta_)
{
  alpha=-alpha_;  delta=delta_;   
}                                 
void calculation::init_beamline(int beamline_)
{
  beamline=beamline_;
  if(beamline==1)  {mode_eps=3;  beta=pi/3;}
  if(beamline==2)  {mode_eps=3;  beta=pi/3;}    
}
void calculation::init_mode_beta(int mode_beta_)
{
  beamline=0;                             
  mode_beta=mode_beta_;
  mode_eps=(mode_beta||mode_sample)?3:1;
}
void calculation::init_mode_sample(int mode_sample_)
{
  beamline=0;                             
  mode_sample=mode_sample_;
  mode_eps=(mode_beta||mode_sample)?3:1;
}
int calculation::n_orders_eff(void)
{
  if(iter.method==3)  return 1;
  else                return n_orders;
}
void calculation::set_polarization(void)
  {light_polarization_mu(e_tth,e_ffi,alpha,delta, eps);}
void calculation::set_polarization(int pol)
{
  if(mode_eps>1) {eps[0]=eps[1]=eps[2]=0;  eps[pol]=1;}
  else           light_polarization_mu(e_tth,e_ffi,alpha,delta, eps);
}
void calculation::free_results(void)
{
  int i;
  for(i=0; i<n_initial_used; i++)  results[i].free();
  if(n_initial>0) delete [] results;  n_initial=n_initial_used=0;
}
void calculation::free_k(void)
{
  if(n_k>0)
  if(k_flag==1)  {delete [] kk;  delete [] k;}  else
  if(k_flag==2)  delete [] kc;
  n_k=0;
}
void calculation::free_orders(void)
{
  if(n_orders>0) {n_orders=0;  delete [] orders;}
}
void calculation::free(void)
{
  integ=n_f=0;
  close_file(fout);     fout=NULL;
  close_file(wf_fin);   wf_fin=NULL;
  close_file(wf_fout);  wf_fout=NULL;
}
void calculation::alloc(int n_initial_)
{
  free_results();
  n_initial=n_initial_*mode_eps;
  results=new array [n_initial];
  int i;   for(i=0; i<n_initial; i++) results[i].n=0;
}
complex *calculation::wave_function(complex *wf)
{
  int i,j;
  if(select_wf<0) {
    if(select_wf==-1)  return wf=wave_function_tmat(wf);      
    if(n_f>0) for(i=0; i<n_f; i++) wf=wave_function(i,0,wf);  
  }
  if(n_ang>0)  for(i=0; i<final.n_th; i++)     
               for(j=0; j<final.n_fi; j++)     
                 wf=wave_function(i,j,wf);
  if(n_mesh>0) for(i=0; i<final.n_1; i++)      
               for(j=0; j<final.n_2; j++)      
                 wf=wave_function(i,j,wf);
  return wf;
}
int calculation::evaluation(int initial)
{
  int o,no, next;
  complex *wf;
  if(n_orders==0)  on_error("calc::evaluation", "undefined orders of scat.");
  n_atoms_last=coord.n;                       
  polarization_flag=0;
  if(iter.method==-1)  results[initial].alloc(n_final * n_c*(max_order+1));
  else                 results[initial].alloc(n_final * n_c*n_orders_eff());
  n_initial_used++;
 
  wf=results[initial].a;
  if(iter.method==3) {                        
    iter.exact();
    wave_function(wf);
    return 0;
  }
  for(o=0, no=0, next=0; o<=max_order; o++) { 
    if(next) iter.next_scattering();          
    if(wf_fin!=NULL && next==0)               
          next=iter.read(wf_fin, o);
    else  next=1;
    if(wf_fout!=NULL)                         
      iter.write(wf_fout, o, max_order);
    if(orders[no]==o || iter.method==-1) {    
      wf=wave_function(wf);
      if(orders[no]==o)  no++;
  } }
  if(iter.method==-1)  eval_recursion(initial);
  return 0;
}
void calculation::eval_recursion(int initial)
{
  int i,j, o,no, nn=n_final*n_c;
  array *results_temp;
  recurs  *val;
  val=new recurs [nn];
  results_temp=new array [nn];
  for(i=0; i<nn; i++) {
    results_temp[i].n=0;  results_temp[i].alloc(max_order+1);
  }
  for(o=j=0; o<=max_order; o++)
  for(i=0; i<nn; i++, j++)  results_temp[i].a[o]=results[initial].a[j];
  results[initial].free();
  for(i=0; i<nn; i++)  val[i].init(max_order, results_temp[i].a);
  for(i=0; i<nn; i++)  results_temp[i].free();
  delete [] results_temp;
  results[initial].alloc(nn * n_orders_eff());
  for(no=j=0; no<n_orders_eff(); no++)
  for(i=0; i<nn; i++, j++)
    results[initial].a[j]=val[i].eval(1,orders[no],0);
  for(i=0; i<nn; i++)  val[i].free();
  delete [] val;
}
void calculation::eval_polarization(void)
{
  int i,j,no,c,index,initial,initial3;  
  numero tt,ff,rr, beamx,beamy,beamz;   
  numero beamxp,beamyp,beamzp;          
  numero beamxi,beamyi,beamzi;          
                                        
  if(mode_eps>1) {                      
    polarization_flag=1;
    for(initial=0; initial<n_initial_used/3; initial++)
    for(no=index=0; no<n_orders_eff(); no++)
    for(i=0; i<final.n_th; i++)
    for(j=0; j<final.n_fi; j++) {initial3=3*initial;
      tt=final.th_out[i];  ff=final.fi[j];
      if(beamline==0) {                 
        if(mode_sample) {
          if(mode_beta)  {e_th=tt-beta;  e_fi=ff;}
          else           {e_th=-e_tth;   e_fi=ff;}
          normalize_angles(e_th, e_fi);
          light_polarization_mu(e_th, e_fi, alpha, delta, eps);
        }  else  {e_th=e_tth;  e_fi=e_ffi;}
      }  else
      if(beamline==1) {                 
        cartesian_to_spherical(sin(beta)*sin(tt)*cos(ff)-cos(beta)*sin(ff),
                               sin(beta)*sin(tt)*sin(ff)+cos(beta)*cos(ff),
                               sin(beta)*cos(tt),
                               rr,e_th,e_fi);    
        e_th=pi/2-e_th;                          
        e_fi=pi  +e_fi;                          
        light_polarization_mu(e_th,e_fi,0,0, eps);
        
      }
      if(beamline==2) {                 
        beamx=-sin(beamline_th)*cos(beamline_fi);
        beamy=-sin(beamline_th)*sin(beamline_fi);
        beamz= cos(beamline_th);
        beamxp= beamx*cos(tt)+beamz*sin(tt);
        beamyp= beamy;
        beamzp=-beamx*sin(tt)+beamz*cos(tt);
        beamxi= beamxp*cos(ff)+beamyp*sin(ff);
        beamyi=-beamxp*sin(ff)+beamyp*cos(ff);
        beamzi= beamzp;
        cartesian_to_spherical(beamxi,beamyi,beamzi,
                               rr,e_th,e_fi);    
        normalize_angles(e_th, e_fi);
        light_polarization_mu(e_th,e_fi,alpha,delta, eps);
      }
      for(c=0; c<n_c; c++, index++)
        results[initial].a[index]=  eps[0]*results[initial3  ].a[index]
                                  + eps[1]*results[initial3+1].a[index]
                                  + eps[2]*results[initial3+2].a[index];
} } }
void calculation::intensity(numero *intensity_,
                            int no, int i, int j, complex *wf)
{
  int c, initial, n_initial_, index;
  complex *wf_;
  if(select_wf<0)     index=(no*n_final+n_ang+n_mesh       +i)*n_c;  else
  if(select_wf<1000)  index=(no*n_final+       i*final.n_fi+j)*n_c;  else
                      index=(no*n_final+n_ang+ i*final.n_2 +j)*n_c;
  if(polarization_flag && select_wf==0)  n_initial_=n_initial_used/3;  
  else                                   n_initial_=n_initial_used;
  if(n_intensities==0 && n_initial_>1)
    on_error("calculation::intensity", "amplitude for multi-initial state");
  if(n_intensities==0) {                                  
    wf_=results[0].a + index;
    for(c=0; c<n_c; c++) wf[c]=wf_[c];
  }  else  {
    for(c=0; c<n_intensities; c++) intensity_[c]=0;
    for(initial=0; initial<n_initial_; initial++) {       
      wf_=results[initial].a + index;
      intensity(intensity_, wf_);        
  } }
  if(select_wf==0 && n_intensities==1)                    
    intensity_[0]=final.transmission[i] * intensity_[0];  
}
void calculation::intensity(numero *II_, complex *wf)
{
  switch(select_wf) {
    case    0: II_[0]+= electron.intensity(wf);  break;  
    case 1000: II_[0]+= electron.intensity(wf);  break;  
    case   -2: II_[0]+= sqr(pi/electron.area)
                        /(sqr(real(scat.k))/2-scat.V0)
                        /(cos(electron.leed_th)*cos(0))  
                        *electron.intensity(wf); break;  
} }
complex *calculation::wave_function(int i, int j, complex *wf)
{
  switch(select_wf) {
    case    0:  electron.wf_far(i,j,wf);   break;   
    case 1000:  electron.wf_near(i,j,wf);  break;   
    case   -2:  electron.wf_LEED(wf);      break;   
  }
  return wf+n_c;                      
}
void calculation::write(int ik)                           
{
  eval_polarization();
  if(n_ang>0) {
    if(verbose)    write_ang(foutput,ik);
    if(fout!=NULL) write_ang(fout,ik);
  } else
  if(n_mesh>0) {
    if(verbose)    write_mesh(foutput,ik);
    if(fout!=NULL) write_mesh(fout,ik);
  } else
  if(n_f>0) {
    if(verbose)    write_f(foutput,ik);
    if(fout!=NULL) write_f(fout,ik);
  }
  free_results();
}
void calculation::write_numeros(FILE *fout_, numero val)
{
  if(ABS(val)<100 && ABS(val)>0.01)  fprintf(fout_," %12.8f", val);
  else                               fprintf(fout_," %12g", val);
}
void calculation::write_numeros(FILE *fout_, numero *intensity_, complex *wf)
{
  int i,c;
  if(n_intensities==0)    
    for(c=0; c<n_c; c++) {
      fprintf(fout_, "  ");
      write_numeros(fout_, real(wf[c]));
      write_numeros(fout_, imag(wf[c]));
    }
  else  for(i=0; i<n_intensities; i++)           
          write_numeros(fout_, intensity_[i]);
}
numero calculation::IIIthfi(int no, numero th, numero fi)
{
  numero ii,ii0, xth,xfi;
  numero th0=final.th[0], th1=final.th[final.n_th-1];
  numero fi0=final.fi[0], fi1=final.fi[final.n_fi-1];
  int ith,ifi,ij;
  while(fi<0)    fi+=2*pi;
  while(fi>2*pi) fi-=2*pi;
  xth=(final.n_th-1)*(th-th0)/(th1-th0);   ith=int(floor(xth-0.001));
  xfi=(final.n_fi-1)*(fi-fi0)/(fi1-fi0);   ifi=int(floor(xfi-0.001));
  if(ifi==-1)  ifi=0;
  if(0<=ith && ith<final.n_th-1 && 0<=ifi && ifi<final.n_fi-1) {
    ij=no*n_ang+ith*final.n_fi+ifi;
    ii0=III0[ij];
    ii=ii0 + (xth-ith)*(III0[ij+final.n_fi]-ii0) + (xfi-ifi)*(III0[ij+1]-ii0);
  }  else  ii=0;
  if(ii<0)  ii=0;
  return ii;
}
numero calculation::IIIave(int no, numero th, numero fi)
{
  if(thave<=1e-6)  return IIIthfi(no,th,fi);
  int  i,j, nn=10, mm=50;
  numero  tth,ffi,val=0, r,f;
  for(i=0; i<nn; i++)
  for(j=0; j<mm; j++)  {
    r=i*thave/nn;
    f=j*2*pi/mm;
    tth=th+r*cos(f);
    ffi=fi+r*sin(f)/cos(pi*th/180);
    val+=IIIthfi(no,tth,ffi);
  }
  return val/(nn*mm);
}
void calculation::write_ang(FILE *fout_, int ik)          
{
  int no,nno,i,j;        
  numero h,fc;           
  complex *wf;           
  numero intensity_[1], val, *III;
  wf=new complex [n_c];
  if(integ==0) {
    nno=n_orders_eff();
    if(thave>1e-6 && n_intensities>0 && final.n_fi>1 && final.n_th>1) {
      if(n_intensities>1)
        on_error(foutput, "scan pd/ms", "only scalars can be averaged");
      III=new numero [nno*n_ang];
      III0=new numero [nno*n_ang];
      for(i=0; i<final.n_th; i++)      
      for(j=0; j<final.n_fi; j++)
      for(no=0; no<nno; no++) {
          intensity(intensity_, no,i,j, wf);
          III0[no*n_ang+i*final.n_fi+j]=intensity_[0];
          III[no*n_ang+i*final.n_fi+j]=0;
      }
      for(no=0; no<nno; no++)          
      for(i=0; i<final.n_th; i++)
      for(j=0; j<final.n_fi; j++)
        III[no*n_ang+i*final.n_fi+j]=IIIave(no, final.th[i], final.fi[j]);
      delete [] III0;
    }
    for(i=0; i<final.n_th; i++)
    for(j=0; j<final.n_fi; j++) {
      write_k(fout_,ik);
      if(final.n_th>1)   fprintf(fout_, " %6.2f", final.th_out[i]*180/pi);
      if(final.n_fi>1)   fprintf(fout_, " %7.2f", final.fi[j]*180/pi);
      for(no=0; no<n_orders_eff(); no++) {
        if(thave>1e-6 && n_intensities>0 && final.n_fi>1 && final.n_th>1)
          write_numeros(fout_, III[no*n_ang+i*final.n_fi+j]);
	else {
          intensity(intensity_, no,i,j, wf);
          write_numeros(fout_, intensity_, wf);
      } }
      fprintf(fout_, "\n");
    }
    if(thave>1e-6 && n_intensities>0 && final.n_fi>1 && final.n_th>1)
      delete [] III;
  } else
  if(integ==1) {
    h = ABS((final.fi[final.n_fi-1]-final.fi[0])/(final.n_fi-1));
    for(i=0; i<final.n_th; i++) {
      write_k(fout_,ik);
      if(final.n_th>1)  fprintf(fout_, " %6.2f", final.th_out[i]*180/pi);
      for(no=0; no<n_orders_eff(); no++) {
        val=0;
        for(j=0; j<final.n_fi; j++) {
          if((j==0 || j==final.n_fi-1) && final.n_fi>1)
                {intensity(intensity_, no,i,j, wf);  val+=0.5*intensity_[0];}
          else  {intensity(intensity_, no,i,j, wf);  val+=    intensity_[0];}
        }
        intensity_[0]=val*h;  write_numeros(fout_, intensity_, wf);
      }
      fprintf(fout_, "\n");
  } } else
  if(integ==2) {
    write_k(fout_,ik);
    h = ABS( (final.th_out[final.n_th-1]-final.th_out[0])/(final.n_th-1)
            *(final.fi[final.n_fi-1]-final.fi[0])/(final.n_fi-1));
    for(no=0; no<n_orders_eff(); no++) {
      val=0;
      for(i=0; i<final.n_th; i++)
      for(j=0; j<final.n_fi; j++) {
        if((i==0 || i==final.n_th-1) && final.n_th>1) fc=0.5;  else  fc=1;
        if((j==0 || j==final.n_fi-1) && final.n_fi>1) fc=0.5*fc;
        intensity(intensity_, no,i,j, wf);
        val+=fc*sin(final.th_out[i])*intensity_[0];
      }
      intensity_[0]=val*h;  write_numeros(fout_, intensity_, wf);
    }
    fprintf(fout_, "\n");
  }
  delete [] wf;
}
void calculation::write_mesh(FILE *fout_, int ik)         
{
  complex *wf;  wf=new complex [n_c];
  numero r,theta,phi, intensity_[1];
  int no,i,j,ij;
  for(i=ij=0; i<final.n_1; i++)
  for(j=0; j<final.n_2; j++, ij++) {
    write_k(fout_,ik);
    if(final.mesh_flag==1) {
      if(final.n_1>1)  fprintf(fout_, " %9.4f", final.x[ij]*a0_au);
      if(final.n_2>1)  fprintf(fout_, " %9.4f", final.z[ij]*a0_au);
    }
    if(final.mesh_flag==2) {
      cartesian_to_spherical(final.x[ij],final.y[ij],final.z[ij],
                             r,theta,phi);
      if(final.n_1>1)  fprintf(fout_, " %6.2f", theta*180/pi);
      if(final.n_2>1)  fprintf(fout_, " %7.2f", phi*180/pi);
    }
    for(no=0; no<n_orders_eff(); no++) {
      intensity(intensity_, no,i,j, wf);
      write_numeros(fout_, intensity_, wf);
    }
    fprintf(fout_, "\n");
  }
  delete [] wf;
}
void calculation::write_f(FILE *fout_, int ik)            
{
  complex *wf;  wf=new complex [n_c];
  numero *intensity_;
  int no,i;
  if(n_intensities>0)  intensity_=new numero [n_intensities];
  for(i=0; i<n_f; i++) {                                  
    write_k(fout_,ik);
    for(no=0; no<n_orders_eff(); no++) {
      intensity(intensity_, no,i,0, wf);
      write_numeros(fout_, intensity_, wf);
    }
    fprintf(fout_, "\n");
  }
  delete [] wf;
}
void calculation::write_k(FILE *fout_, int ik)            
{
  if(n_k>1 || ms)
  if(k_flag==1 || ms)  fprintf(fout_, " %9.4f", kk[ik]);  else
  if(k_flag==2)  fprintf(fout_, " %9.4f %9.4f",
                                real(kc[ik])/units_conv_fct(k_units),
                                imag(kc[ik])/units_conv_fct(k_units));
}
void calculation::write_ang(FILE *fout_, char *title)     
{
  write_title(fout_, title);
  if(final.n_th>1 && (integ==0 || integ==1))  fprintf(fout_,"  theta");
  if(final.n_fi>1 && integ==0)                fprintf(fout_,"     phi");
  write_orders(fout_);
}
void calculation::write_mesh(FILE *fout_, char *title)    
{
  write_title(fout_, title);
  if(final.mesh_flag==1) {
    if(final.n_1>1)  fprintf(fout_,"  x1 (Ang)");
    if(final.n_2>1)  fprintf(fout_,"  x2 (Ang)");
  }
  if(final.mesh_flag==2) {
    if(final.n_1>1)  fprintf(fout_,"  theta");
    if(final.n_2>1)  fprintf(fout_,"     phi");
  }
  write_orders(fout_);
}
void calculation::write_f(FILE *fout_, char *title)       
{
  write_title(fout_, title);
  write_orders(fout_);
}
void calculation::write_title(FILE *fout_, char *title)   
{
  fprintf(fout_, "--- %s\n", title);
  if(n_k>1 || ms)
  if(k_flag==1 || ms)  fprintf(fout_,"     %5.5s", str_units(k_units));  else
  if(k_flag==2)  fprintf(fout_,"  complex k in %5.5s", str_units(k_units));
}
void calculation::write_orders(FILE *fout_)               
{
  int i, c;
  if(iter.method==3)  fprintf(fout_,"        exact");     
  else
  for(i=0; i<n_orders_eff(); i++)
    if(n_intensities==0)
      for(c=0; c<n_c; c++)
        fprintf(fout_,"      order %2d   component %d", orders[i], c+1);
    else
      for(c=0; c<n_intensities; c++)
        if(i==0 && c==0)  fprintf(fout_,"     order %2d", orders[i]);
        else  if(c==0)    fprintf(fout_," --- order %2d", orders[i]);
              else        fprintf(fout_,"-------------", orders[i]);
  fprintf(fout_,"\n");
}
void calculation::init_wf_tmat(int c, int l, int m)
{
  int a;
  vector gg;
  scat.translation_type=1;           
  for(a=0; a<coord.n; a++) if(coord.atom[a]>=0) {
    gg.alloc(l,n_c);
    gg.re[c*sqr(l+1)+l*l+l+m]=1;     
    if(a!=center) {
      scat.get_ij(a,center);         
      scat.get_ij();                 
      scat.propagate(gg,lmax_inter);
    }
    scat.R_1TR(gg, a);               
    iter.o.g[a].copy(gg);
    iter.n.g[a].copy(gg);
    gg.free();
  }  else  {iter.o.g[a].alloc(-1,n_c);  iter.n.g[a].alloc(-1,n_c);}
  scat.translation_type=0;           
}
complex *calculation::wave_function_tmat(complex *wf)
{
  int a;
  vector gg, sum;
  scat.translation_type=1;           
                                     
  sum.alloc(-1,n_c);                 
  for(a=0; a<coord.n; a++) if(coord.atom[a]>=0) {
    gg.copy(iter.n.g[a]);
    if(a!=center) {
      scat.get_ij(center,a);         
      scat.get_ij();                 
      scat.propagate(gg,lmax_eff);   
    }                                
    sum.add(gg);                     
    gg.free();
  }
  scat.translation_type=0;           
  int c,i,lm,ll=lp12[sum.lmax];
  for(c=i=0; c<n_c ; c++)            
  for(lm=0; lm<lp12[lmax]; lm++, i++)
    if(lm<ll)  wf[i]=complex(sum.re[c*ll+lm], sum.im[c*ll+lm]);
    else       wf[i]=0;
  sum.free();
  return wf+dim;
}
void calculation::scan_tmat(char *name, int lmax_)
{
  int ik, i,j,ij, c,l,m, cp,lp,mp;   
  numero  val;                       
  complex *aa;                       
  lmax=lmax_;  n_c=scat.get_n_c();
  dim=n_c*sqr(lmax+1);               
  scat.set_lmax();                   
  lmax_inter=lmax_eff;               
  if(lmax>lmax_eff) {                
    free_lmax();  lmax_eff=lmax;
    init_lmax();
  }
  coord.sort(-1);                    
  
  if(coord.get_atom(coord.x0,coord.y0,coord.z0)<0)
    coord.add_atom(coord.x0,coord.y0,coord.z0,0,0,-1);
  center=coord.get_atom(coord.x0,coord.y0,coord.z0);
  
  init_final_f(dim/n_c);  select_wf=-1;  
                                         
  init(name);                            
  iter.alloc();
  
  if(fout!=NULL) {
    fprintf(fout, "%d", n_k);
    if(n_k>1)  fprintf(fout, " E(eV)");
    fprintf(fout, " %d full %d\n", lmax, n_c);
  }
  if(verbose) {
    fprintf(foutput,"--- t matrix of the cluster,  lmax=%d", lmax);
    fprintf(foutput,"  %d component(s),  output file %s\n", n_c, name);
  }
  
  aa=new complex [dim*dim];          
  for(ik=0; ik<n_k; ik++) {          
    val=scat.iimfp.val(k[ik])/2;   
    scat.k=complex(k[ik], val);
    if(iter.method==3) {             
      scat.get_exact(n_c,lmax_inter);
      if(verbose)
        fprintf(foutput,
                "    exact inversion completed with lmax=%d\n", lmax_inter);
    }
    for(cp=j=0; cp<n_c ; cp++)       
    for(lp=0; lp<=lmax; lp++) {
      for(mp=-lp; mp<=lp; mp++, j++) {
        alloc(1);                    
                                     
        if(wf_fin==NULL)
          init_wf_tmat(cp, lp, mp);  
        evaluation(0);               
        iter.clear();                
        for(i=0; i<dim; i++) aa[i*dim+j]=results[0].a[i];
        free_results();
      }
      if(verbose && iter.method!=3)
        fprintf(foutput,"    ik=%2d   c'=%d   l'=%2d\n", ik+1, cp+1, lp);
    }
    if(iter.method==3)  scat.free_exact();     
    scat.free_Green_functions();     
    scat.free_t_matrix();            
  
    if(n_k>1)  fprintf(fout, "\n%g\n\n", sqr(k[ik])/2*au_eV);
    for(c=ij=0; c<n_c ; c++)         
    for(l=0; l<=lmax; l++)
    for(m=-l; m<=l; m++)
    for(cp=0; cp<n_c ; cp++)
    for(lp=0; lp<=lmax; lp++)
    for(mp=-lp; mp<=lp; mp++, ij++) {
      if(n_c>1)  fprintf(fout, "%d",        c+1);
                 fprintf(fout, " %d %2d",   l,m);
      if(n_c>1)  fprintf(fout, " %d",       cp+1);
                 fprintf(fout, " %d %2d",   lp,mp);
      if(mod(aa[ij])<1e-20)  fprintf(fout, " 0 0\n");  else
                 fprintf(fout, " %g %g\n",  real(aa[ij]), imag(aa[ij]));
  } }
  if(coord.atom[center]<0)           
    coord.delete_atom(coord.x0,coord.y0,coord.z0);
  free();                            
  iter.free();                       
  scat.free_rotations();             
  free_lmax();                       
  delete [] aa;
}
green_z::green_z(void)
{
  re = new numero [size_of_green_z];
  im = new numero [size_of_green_z];
}
void green_z::init(complex k, numero dd, int translation_type, int att_type)
{
  int i;  int m1,l1,l2,l3, il;
  numero fct=1,j,y,hr,hi, *pre,*pim, *spherical_hr, *spherical_hi, kdr,kdi;
  complex kd, val;
  kd=k*dd;  d=dd;  trans=translation_type;
  kdr=real(kd);  kdi=imag(kd);
  spherical_hr=new numero [2*lmax_eff+1];
  spherical_hi=new numero [2*lmax_eff+1];
  if(att_type==0)  for(l3=0; l3<=2*lmax_eff; l3++) {
                     val=besselilh(l3,kd);
                     spherical_hr[l3]=real(val);   
                     spherical_hi[l3]=imag(val);
                   }
  else {
    fct = exp(-kdi);                            
                                                
    if(kdr)  fct *= sqrt(pi/(2*kdr));           
                                                
    for(l3=0; l3<=2*lmax_eff; l3++) if(kdr) {   
      besselJY(l3+0.5,kdr,j,y,hr,hi);  il=l3%4; 
      if(translation_type==1) {y=-j; j=0;}      
      if     (il==0) {hr=-y; hi= j;}      
      else if(il==1) {hr=-j; hi=-y;}      
      else if(il==2) {hr= y; hi=-j;}      
      else           {hr= j; hi= y;}      
      spherical_hr[l3]=fct*hr;            
      spherical_hi[l3]=fct*hi;
    } else {spherical_hr[l3]=1; spherical_hi[l3]=0;}
  }
  for(l1=0, pre=re,pim=im, i=0; l1<=lmax_eff; l1++)
  for(m1=-l1; m1<=l1; m1++)
  for(l2=(m1>0?m1:-m1); l2<=lmax_eff; l2++, pre++,pim++)
  for(l3=(l1>l2?l1-l2:-l1+l2), *pre=*pim=0; l3<=l1+l2; l3+=2) {
    (*pre)+=spherical_hr[l3]*Gaunt_integrals[i];   
    (*pim)+=spherical_hi[l3]*Gaunt_integrals[i++]; 
  }                                                
                                                   
  delete [] spherical_hr;
  delete [] spherical_hi;
}                        
int vector::Gz(green_z *G, int lmax_final)
{
  int n,c,nele,nelet, l1,l2,m1, *llpl;  vector temp;
  numero *Gr,*Gi,*ptr,*pti, *rre,*iim;
  if(lmax<0) return 0;                
  if(lmax_final>lmax_eff)  lmax_final=lmax_eff;
  temp.alloc(lmax_final, n_c);
  llpl=new int [lmax+1];   for(l2=0; l2<=lmax; l2++)  llpl[l2]=l2*l2+l2;
  for(c=0; c<n_c; c++) {
    nele =c*lp12[lmax];
    nelet=c*lp12[temp.lmax];
    ptr=temp.re+nelet;   rre=re+nele;
    pti=temp.im+nelet;   iim=im+nele;
    for(l1=0, Gr=G->re,Gi=G->im; l1<=temp.lmax; l1++)
    for(m1=-l1; m1<=l1; m1++, ptr++,pti++) { *ptr=*pti=0;
      for(l2=(m1>0?m1:-m1); l2<=lmax; l2++, Gi++,Gr++) {
        n=llpl[l2]+m1;
        (*ptr)+= (*Gr)*rre[n]-(*Gi)*iim[n];
        (*pti)+= (*Gr)*iim[n]+(*Gi)*rre[n];
      }
      if(l2<=lmax_eff) {Gi+=(lmax_eff-l2+1); Gr+=(lmax_eff-l2+1);}
  }}
  free();  delete [] llpl;
  *this=temp;
  return 0;
}
rotation_th::rotation_th(void)  {val = new numero [size_of_rotation_th];}
int rotation_th::init(numero tth)
{
  int l,m1,m2,t,e1,e2;   numero qsi,eta,pot, *pval,*ppval,*pnval;
  qsi=cos(tth/2.0);  eta=sin(tth/2.0);    
                                          
  for(l=0, pval=val; l<=lmax_eff; l++) {  
    pnval=pval;                           
    for(m1=-l;  m1<=l; m1++)              
    for(m2=0, pval+=l; m2<=l; m2++, pval++)
    for(t=(0>m1-m2?0:m1-m2), *pval=0; t<=(l+m1<l-m2?l+m1:l-m2); t++) {
      e2=t+t-m1+m2;  e1=l+l-e2;
      if(!qsi) if(e1) pot=0; else pot=1; else pot=pow(qsi,e1+0.0);
      if(!eta) if(e2) pot=0; else pot=pot; else pot=pot*pow(eta,e2+0.0);
      *pval+= (t&1?-1:1) * pot
              * (fact2[m1+l]/fact[l+m1-t]) * (fact2[l-m1]/fact[l-m2-t])
              * (fact2[m2+l]/fact[t])      * (fact2[l-m2]/fact[t-m1+m2]);
    }
    for(m1=-l, ppval=pval-1;  m1<=l; m1++) {        
      for(m2=-l;  m2<0; m2++, pnval++,ppval--)      
        *pnval=((m1&1)^(m2&1))?-(*ppval):(*ppval);  
        pnval+=(l+1);  ppval-=(l+1);                
    }                                               
  }                                                 
                                                    
  return 0;                                         
}
int vector::Rth(rotation_th *R, int op)
{
  int l,m1,m2, c,nele,nelet;  vector temp;
  numero *rlm,*pre,*pim,*ptre,*ptim,*tr,*ti,ee;
  if(lmax<0) return 0;            
  if(*(R->val)>1) {
    if(*(R->val)==2) return 0;    
    for(c=0; c<n_c; c++)          
    for(l=0, pre=re+(c*lp12[lmax]),
             pim=im+(c*lp12[lmax]); l<=lmax; l++) {  
      pre=pre+l; pim=pim+l; tr=pre; ti=pim;          
      for(m1=0; m1<=l; m1++, pre++,pim++, tr--,ti--)
      if((l&1)^(m1&1)) {ee=*pre;  *pre=-(*tr);  (*tr)=-ee;
                        ee=*pim;  *pim=-(*ti);  (*ti)=-ee;}
      else             {ee=*pre;  *pre=  *tr;   (*tr)= ee;
                        ee=*pim;  *pim=  *ti;   (*ti)= ee;}
    }
    return 0;
  }
  temp.alloc(lmax,n_c);
  for(c=0; c<n_c; c++) {          
    nele =c*lp12[lmax];
    nelet=c*lp12[temp.lmax];
    ptre=temp.re+nelet;   pre=re+nele;
    ptim=temp.im+nelet;   pim=im+nele;  rlm=R->val;
    if(op)
    for(l=0; l<=temp.lmax; l++) { tr=ptre;  ti=ptim;
    for(m1=-l; m1<=l; m1++, pre++,pim++) { ptre=tr;  ptim=ti;
    for(m2=-l; m2<=l; m2++, rlm++, ptre++,ptim++)  {
      (*ptre)+=(*rlm)*(*pre);
      (*ptim)+=(*rlm)*(*pim);
    }}}
    else
    for(l=0; l<=temp.lmax; l++) { tr=pre;  ti=pim;
    for(m1=-l; m1<=l; m1++, ptre++,ptim++) { pre=tr;  pim=ti;
    for(m2=-l; m2<=l; m2++, rlm++, pre++,pim++)  {
      (*ptre)+=(*rlm)*(*pre);
      (*ptim)+=(*rlm)*(*pim);
    }}}}
  free();
  *this=temp;
  return 0;
}
rotation_fi::rotation_fi(void)
{
  re = new numero [size_of_rotation_fi];
  im = new numero [size_of_rotation_fi];
}
int rotation_fi::init(numero ffi)
{
  int m;  numero ang;
  if(ABS(ffi)   <coord.eps_th)    {re[0]=2; return 1;}
  if(ABS(ffi-pi)<coord.eps_th ||
     ABS(ffi+pi)<coord.eps_th)    {re[0]=3; return 2;}
  for(m=0, ang=-lmax_eff*ffi; m<=2*lmax_eff; m++) {
    re[m]= cos(ang);                  
    im[m]=-sin(ang);  ang+=ffi;       
  }
  return 0;
}
int vector::Rfi(rotation_fi *R, int op)
{
  int l,m,c;  numero *pre,*pim,*pRre,*pRim,tre,tim;
  if(lmax<0) return 0;           
  if(*(R->re)>1) {
    if((*(R->re)==3||op>1)&&(*(R->re)==2||op<=1)) 
    for(c=0, pre=re, pim=im; c<n_c; c++)          
    for(l=0; l<=lmax; l++)                        
    for(m=-l; m<=l; m++, pre++,pim++)
      if(m&1) {*pre=-(*pre); *pim=-(*pim);}       
    return 0;                                     
  }                                               
                                                  
  for(c=0, pre=re, pim=im; c<n_c; c++)
  for(l=0; l<=lmax; l++)
  for(m=-l, pRre=R->re+(lmax_eff-l), pRim=R->im+(lmax_eff-l); m<=l;
    m++, pRre++,pRim++,pre++,pim++) {
    switch(op) {
       case 0:  case 2: tre=(*pRre)*(*pre)-(*pRim)*(*pim);
       	       	        tim=(*pRre)*(*pim)+(*pRim)*(*pre);  break;
       case 1:  case 3: tre=(*pRre)*(*pre)+(*pRim)*(*pim);
                        tim=(*pRre)*(*pim)-(*pRim)*(*pre);  break;
    }
    if(op>1&&(m&1))  {*pre=-tre; *pim=-tim;}
    else             {*pre= tre; *pim= tim;}
  }
  return 0;
}
vector::vector(void)  {lmax=-1;  n_c=1;}
void vector::free(void)                          
{
  if(lmax>=0) {delete [] re;  delete [] im;}
  lmax=-1;  n_c=1;
}
void vector::alloc(int lmax_, int n_c_)          
{
  int i,dimension;
  if(lmax>=0)  free();
  n_c = n_c_;                      
  lmax=lmax_;                      
  if(lmax>lmax_eff) lmax=lmax_eff;
  if(lmax>=0) {
    dimension=n_c*lp12[lmax];      
    re = new numero [dimension];
    im = new numero [dimension];
    for(i=0; i<dimension; i++)  re[i]=im[i]=0;
} }
void vector::copy_from_array(complex *a, int lmax_, int n_c_)
{
  int n,c;  numero *rre,*iim;  complex *va;
  free();  alloc(lmax_,n_c_);  rre=re; iim=im;  va=a;
  if(lmax>=0)
  for(c=0; c<n_c; c++)
  for(n=0; n<lp12[lmax]; n++, rre++, iim++, va++) {
    *rre=real(*va);  *iim=imag(*va);
} }
void vector::copy_to_array(complex *a, int lmax_, int cc)
{
  int l,m,c;  numero *rre,*iim;  complex *va;
  rre=re; iim=im;  va=a;
  if(lmax>=0)
  for(c=0; c<n_c; c++)
  for(l=0; l<=lmax; l++)
  for(m=-l; m<=l; m++, rre++, iim++)  if((c==cc || cc<0) && l<=lmax_)  {
    *va=complex(*rre, *iim);
    va++;
} }
void vector::copy_to_array(complex *a, int cc)  {copy_to_array(a, lmax, cc);}
void vector::copy_to_array(complex *a)  {copy_to_array(a,-1);}
void vector::copy(vector &a, int lmax_, int n_c_, int c_flag)
{
  int n,c;  numero *tr,*ti,*rre,*iim;
  n_c_=(a.n_c<n_c_)?n_c_:a.n_c;
  lmax_=(a.lmax>lmax_)?lmax_:a.lmax;
  free();  alloc(lmax_,n_c_);
  if(lmax>=0)
  for(c=0, tr=re, ti=im; c<n_c_; c++) {
    if((c%a.n_c)==0) {rre=a.re; iim=a.im;}
    for(n=0; n<lp12[a.lmax]; n++, rre++, iim++)
      if(n<lp12[lmax]) {
        if(c_flag<0 || c==c_flag)  {*tr=*rre; *ti=*iim;}
        else                       {*tr=0;    *ti=0;}
        tr++; ti++;
} }   }
void vector::copy(vector &a)                       
{
  copy(a, lmax_eff, a.n_c, -1);
}
void vector::substitute_lmax(int lmax_)            
{
  vector temp;  temp.copy(*this, lmax_, n_c, -1);   free();  *this=temp;
}
void vector::substitute_n_c(int n_c_, int c_flag)  
{
  vector temp;  temp.copy(*this, lmax, n_c_, c_flag);  free();  *this=temp;
}
int vector::add(vector &b)
{
  int l,m, c;  numero *ar,*ai,*br,*bi,*tr,*ti;  vector temp;
  if(b.lmax<0) return 0;
  if(lmax<0) {copy(b);  return 0;}
  if(n_c!=b.n_c)  on_error(foutput,"vector::add","unphysical case");
  temp.alloc((b.lmax>lmax)?b.lmax:lmax,n_c);
  for(c=0, tr=temp.re, ti=temp.im,
      ar=re, ai=im, br=b.re, bi=b.im; c<temp.n_c; c++) {
    for(l=0; l<=temp.lmax; l++)
    for(m=-l; m<=l; m++, tr++, ti++) {
      (*tr)=(*ti)=0;
      if(l<=lmax) {(*tr)+=(*ar); (*ti)+=(*ai); ar++; ai++;}
      if(l<=b.lmax) {(*tr)+=(*br); (*ti)+=(*bi); br++; bi++;}
  } }
  free();  *this=temp;  return 0;
}
int vector::ext(numero x)
{
  int l,m, c;  numero *br,*bi;
  if(lmax<0 || x==1) return 0;
  for(c=0, br=re, bi=im; c<n_c; c++) {
    if(x==-1) for(l=0; l<=lmax; l++)
              for(m=-l; m<=l; m++, br++,bi++) {(*br)=-(*br); (*bi)=-(*bi);}
    else      for(l=0; l<=lmax; l++)
              for(m=-l; m<=l; m++, br++,bi++) {(*br)=x*(*br); (*bi)=x*(*bi);}
  }
  return 0;
}
int vector::tmat(t_matrix *t)            
{
  if(lmax<0)  return 0;
  if(t->ta<0)  {free();  return 0;}
  if(lmax>t->lmax)  substitute_lmax(t->lmax);
  if(t->select==0)  tmat_diagonal(t);   else
  if(t->select==1)  tmat_full(t);       else
  if(t->select==2)  tmat_full(t);       else
  if(t->select==4)  tmat_so(t);         else
                    on_error(foutput,"tmat", "undefined case");
  return 0;
}
void vector::tmat_diagonal(t_matrix *t)  
{                                        
  int l,m,c,j;                           
  numero *br,*bi,*tr,*ti,rr,ii,ttr,tti;  
                                         
  if(t->n_c!=1 && n_c!=1 && t->n_c!=n_c) 
    on_error(foutput,"vector::tmat_diagonal","unphysical case");
  for(c=0, br=re, bi=im; c<n_c; c++) {
    if(t->n_c==1 || c==0)  {tr=t->re; ti=t->im;}
    else {
      tr=t->re+c*(t->lmax+1);
      ti=t->im+c*(t->lmax+1);
    }
    for(l=0; l<=lmax; l++, tr++, ti++) { 
      ttr=*tr;  tti=*ti;                 
      if(t->n_c>1 && n_c==1) {
        for(j=0; j<t->n_c; j++) {ttr+=(*(tr+j*(t->lmax+1)));
                                 tti+=(*(ti+j*(t->lmax+1)));}
        ttr/=t->n_c;   tti/=t->n_c;
      }
      for(m=-l; m<=l; m++, br++,bi++) {
        rr =ttr*(*br);  ii =ttr*(*bi);
        rr-=tti*(*bi);  ii+=tti*(*br);
        *br=rr;  *bi=ii;   
} } } }
void vector::tmat_full(t_matrix *t)      
{                                        
  int c,l,m,cp,lp,mp;
  numero *ar,*ai,*br,*bi,*tr,*ti;
  vector a;                              
                                         
  if(t->n_c!=n_c)
    on_error(foutput,"vector::tmat_full","unphysical case");
  a.alloc(t->lmax,n_c);                  
  ar=a.re;  ai=a.im;  tr=t->re;  ti=t->im;
  for(c=0; c<n_c; c++)
  for(l=0; l<=t->lmax; l++)
  for(m=-l; m<=l; m++) {
    br=re;  bi=im;  *ar=*ai=0;
    for(cp=0; cp<n_c; cp++)
    for(lp=0; lp<=t->lmax; lp++)
    for(mp=-lp; mp<=lp; mp++) {
      if(t->select==2 || m==mp) {        
        if(lp<=lmax) {                   
          *ar+=(*tr)*(*br)-(*ti)*(*bi);
          *ai+=(*tr)*(*bi)+(*ti)*(*br);
        }
        tr++;  ti++;
      }
      if(lp<=lmax)  {br++;  bi++;}
    }
    ar++;  ai++;
  }
  free();  *this=a;
}
void vector::tmat_so(t_matrix *t)
{
  int l,m, lp12p1=lp12[lmax]+1, lp12p0=lp12[lmax];
  numero *ar, *ai, *br,*bi,*tr,*ti,rr,ii,
         truu,trud,trdu,trdd,tiuu,tiud,tidu,tidd;
  vector a;  a.copy(*this);
  for(l=0, tr=t->re,ti=t->im, br=re,bi=im, ar=a.re,ai=a.im; l<=lmax; l++)
  for(m=-l; m<=l; m++, ar++,ai++, br++,bi++, tr+=4,ti+=4) {
    truu=*tr;   trdd=*(tr+2);   tiuu=*ti;   tidd=*(ti+2);
    if(a.n_c==1) {
      *br=((truu+trdd)*(*ar)-(tiuu+tidd)*(*ai))/2;
      *bi=((truu+trdd)*(*ai)+(tiuu+tidd)*(*ar))/2;
    }  else  {
      trud=*(tr+1);   trdu=*(tr+3);   tiud=*(ti+1);   tidu=*(ti+3);
                rr =truu*(*ar)         -tiuu*(*ai);
                ii =truu*(*ai)         +tiuu*(*ar);
      if(m<l)  {rr+=trud*(*(ar+lp12p1))-tiud*(*(ai+lp12p1));
                ii+=trud*(*(ai+lp12p1))+tiud*(*(ar+lp12p1));}
      *br=rr;  *bi=ii;
                rr =trdd*(*(ar+lp12p0))-tidd*(*(ai+lp12p0));
                ii =trdd*(*(ai+lp12p0))+tidd*(*(ar+lp12p0));
      if(m>-l) {rr+=trdu*(*(ar-1))     -tidu*(*(ai-1));
                ii+=trdu*(*(ai-1))     +tidu*(*(ar-1));}
      *(br+lp12p0)=rr;  *(bi+lp12p0)=ii;
  } }
  a.free();
}
int vector::spin_rotation(numero al, numero be, numero ga)
{
  if(lmax<0 || n_c==1)  return 0;            
  numero cp=cos((al+ga)/2), sp=sin((al+ga)/2), cb=cos(be/2),
         cm=cos((al-ga)/2), sm=sin((al-ga)/2), sb=sin(be/2),
         uur,uui,udr,udi,dur,dui,ddr,ddi, ur,dr,ui,di, *br,*bi;
  int n,dim=lp12[lmax];                   
  uur=cp*cb;  uui=-sp*cb;       udr=-cm*sb;  udi=sm*sb;   
  dur=cm*sb;  dui= sm*sb;       ddr= cp*cb;  ddi=sp*cb;   
  for(n=0, br=re, bi=im; n<lp12[lmax]; n++, br++, bi++) {
    ur=(*br);    dr=(*(br+dim));
    ui=(*bi);    di=(*(bi+dim));                   
    (*br)      =(uur*ur-uui*ui)+(udr*dr-udi*di);   
    (*bi)      =(uur*ui+uui*ur)+(udr*di+udi*dr);   
    (*(br+dim))=(dur*ur-dui*ui)+(ddr*dr-ddi*di);
    (*(bi+dim))=(dur*ui+dui*ur)+(ddr*di+ddi*dr);   
  }
  return 0;
}
cluster::cluster(void)
{
  x0=y0=z0=0;
  oriented=0;
  zs=0;  surface=0;
  Rmax=infinity;
  n=nn=0;
  nnn=infinite_int;
  n_emit=0;
  eps_d=0.00005;
  eps_th=0.003*pi/180;
}
int cluster::set_n(void)
{
  n=(nnn>nn)?nn:nnn;
  return n;
}
void cluster::free(void)
{
  free_coordinates();
  free_emitters();
}
void cluster::free_coordinates(void)
{
  if(nn) {
    delete [] x;  delete [] y;  delete [] z;
    delete [] atom;  delete [] atZ;  delete [] th;  delete [] fi;
    n=nn=0;
} }
void cluster::free_emitters(void)            
{
  if(ABS(n_emit)>0) {n_emit=0;  delete [] emit;  delete [] emit_s;}
}
int cluster::read_atomic_type(FILE *fin)
{
  char name[str_length_comments];  read_name(fin,name);
  int  jj=0;
                                             
  jj=element_properties(name);               
          
  if(jj==0)  jj=read_int(name);
  return jj;
}
void cluster::read(FILE *fin, char *my_file)
{
  FILE *fcoordinates;
  if(!strcmpC(my_file,"inline")) fcoordinates=fin;
  else                           fcoordinates=open_file(foutput,my_file,"r");
  free();
  read(fcoordinates);
  if(strcmpC(my_file,"inline"))  fclose(fcoordinates);
}
void cluster::read(FILE *fcoordinates)
{
  int i,j;
  numero factor2;
  char units2[str_length_units];
  nn=read_int(fcoordinates);  set_n();
  read_name(fcoordinates,units2);
  factor2=units_conversion_factor(units2);
  x=new numero [nn];  y=new numero [nn];  z=new numero [nn];
  atom=new int [nn];  atZ=new int [nn];
  th= new numero [nn];  fi= new numero [nn];
  for(j=0; j<nn; j++) {
    i=read_int(fcoordinates)-1;
    atZ[i]=atom[i]=read_atomic_type(fcoordinates);
    x[i]=factor2*read_numero(fcoordinates);
    y[i]=factor2*read_numero(fcoordinates);
    z[i]=factor2*read_numero(fcoordinates);
    if(oriented)  {
      th[i]=pi/180*read_numero(fcoordinates);
      fi[i]=pi/180*read_numero(fcoordinates);
    }  else {th[i]=fi[i]=0;}
} }
void cluster::write(char *name, char *units)
{
  FILE *fout;
  set_n();
  if(!strcmpC(name,"none"))  fout=foutput;
  else                       fout=open_file(foutput,name,"w");
  write(fout, units);
  close_file(fout);
}
void cluster::write(FILE *fout, char *units)
{
  int i;     numero factor=units_conversion_factor(units);
  fprintf(fout, "%3d %s\n", n, units);
  for(i=0; i<n; i++) {
    fprintf(fout, "%3d %3d %15.9f %15.9f %15.9f", i+1, atom[i],
            x[i]/factor, y[i]/factor, z[i]/factor);
    if(oriented)
      fprintf(fout, "   %9.4f %9.4f", 180/pi*th[i], 180/pi*fi[i]);
    fprintf(fout, "\n");
} }
numero cluster::distance(numero x0_, numero y0_, numero z0_,
                         numero x1, numero y1, numero z1)
{
  numero  ddd=sqrt( (x0_-x1)*(x0_-x1)       
                   +(y0_-y1)*(y0_-y1)       
                   +(z0_-z1)*(z0_-z1));     
  if(surface) ddd=ddd+zs-z1;                
                                            
  return ddd;
}
int cluster::get_atom(numero x_, numero y_, numero z_, numero err)
{                                          
  int i;                                   
                                           
  for(i=0; i<nn; i++)
    if(sqr(x[i]-x_)+sqr(y[i]-y_)+sqr(z[i]-z_) < sqr(err))
    return i;
  return -1;
}
int cluster::get_atom(numero x_, numero y_, numero z_)
  {return get_atom(x_,y_,z_,eps_d);}
void cluster::sort(int emitter)          
{                                        
  int *atom_,*atZ_,*am, i,j,l;           
  numero *x_,*y_,*z_,*th_,*fi_,*d, dis, xe,ye,ze;    
  x_=new numero [nn];                    
  y_=new numero [nn];                    
  z_=new numero [nn];                    
  atom_=new int [nn];                    
  atZ_=new int [nn];                     
  th_=new numero [nn];                   
  fi_=new numero [nn];                   
  d=new numero [nn];
  am=new int [nn];
  if(emitter>=0) {                                 
    xe=x[emitter];  ye=y[emitter];  ze=z[emitter]; 
  }  else  {xe=x0; ye=y0; ze=z0;}                  
  for(i=n=0; i<nn; i++) {                
    d[i]=distance(xe,ye,ze,              
                  x[i],y[i],z[i]);       
    if(d[i]<=Rmax) n++;                  
    am[i]=-1;
  }
                                         
  if(nn<n)  n=nn;                        
  if(nnn<n) n=nnn;                       
  for(i=0; i<nn; i++) {
    if(i>=n)   for(l=0; am[l]>=0; l++);
    else       for(j=0,dis=infinity; j<nn; j++)
                 if(am[j]<0 && d[j]<dis) {dis=d[j]; l=j;}
    am[l]=i;
    x_[i]=x[l];  y_[i]=y[l];  z_[i]=z[l];
    atom_[i]=atom[l];  atZ_[i]=atZ[l];  th_[i]=th[l];  fi_[i]=fi[l];
  }
  for(j=0; j<ABS(n_emit); j++)  emit[j]=am[emit[j]];
  delete [] x;  delete [] y;  delete [] z;
  delete [] atom;  delete [] atZ;  delete [] th;  delete [] fi;
  delete [] d;  delete [] am;
  x=x_;  y=y_;  z=z_;  atom=atom_;  atZ=atZ_;  th=th_;  fi=fi_;
}
int cluster::add_atom(numero x_, numero y_, numero z_,
                      numero th_, numero fi_, int at)      
{
  if(distance(x0,y0,z0, x_,y_,z_)>Rmax)  return 0;
  int i,*atom_,*atZ_;  numero *xx,*yy,*zz,*tt,*ff;
  i=get_atom(x_,y_,z_);
  if(i>=0) {atZ[i]=atom[i]=at;  th[i]=th_;  fi[i]=fi_;  return 0;}
  if(nn) {xx=x;  yy=y;  zz=z;  atom_=atom;  atZ_=atZ;  tt=th;  ff=fi;}
  x=new numero [nn+1];  y=new numero [nn+1];  z=new numero [nn+1];
  atom=new int [nn+1];  atZ=new int [nn+1];
  th=new numero [nn+1];  fi=new numero [nn+1];
  for(i=0; i<nn; i++) {
    x[i]=xx[i];  y[i]=yy[i];  z[i]=zz[i];
    atom[i]=atom_[i];  atZ[i]=atZ_[i];  th[i]=tt[i];  fi[i]=ff[i];
  }
  x[nn]=x_;  y[nn]=y_;  z[nn]=z_;
  atZ[nn]=atom[nn]=at;
  th[nn]=th_;  fi[nn]=fi_;
  if(nn) {
    delete [] xx;  delete [] yy;  delete [] zz;
    delete [] atom_;  delete [] atZ_;  delete [] tt;  delete [] ff;
  }
  nn++;  set_n();  return 0;
}
int cluster::delete_atom(numero x_, numero y_, numero z_)  
{
  if(distance(x0,y0,z0, x_,y_,z_)>Rmax)  return 0;
  int ii=get_atom(x_,y_,z_);  if(ii<0)  return 0;
  if(nn==1) {free();  return 0;}
  int i,j,*atom_,*atZ_;  numero *xx,*yy,*zz,*tt,*ff;
  xx=x;  yy=y;  zz=z;  atom_=atom;  atZ_=atZ;  tt=th;  ff=fi;
  x=new numero [nn-1];  y=new numero [nn-1];  z=new numero [nn-1];
  atom=new int [nn-1];  atZ=new int [nn-1];
  th=new numero [nn-1];  fi=new numero [nn-1];
  for(i=j=0; i<nn; i++) if(i!=ii) {
    x[j]=xx[i];  y[j]=yy[i];  z[j]=zz[i];
    atom[j]=atom_[i];  atZ[j]=atZ_[i];  th[j]=tt[i];  fi[j]=ff[i];
    j++;
  }
  delete [] xx;  delete [] yy;  delete [] zz;
  delete [] atom_;  delete [] atZ_;  delete [] tt;  delete [] ff;
  nn=nn-1;  set_n();  return 0;
}
void cluster::add_layer(numero x_, numero y_, numero z_, numero a, numero b,
                        numero al, numero be, numero th_, numero fi_, int at)
{
  int i,j,imax,jmax;
  numero d, ca=a*cos(al),cb=b*cos(al+be), sa=a*sin(al),sb=b*sin(al+be);
  if((Rmax/infinity)>0.9)  on_error(foutput,"add_layer","Rmax too large");
  d=Rmax + sqrt(sqr(x_-x0)+sqr(y_-y0));
  imax=int(floor(d/a+1));  jmax=int(floor(d/b+1));
  for(i=-imax; i<=imax; i++)
  for(j=-jmax; j<=jmax; j++)
    add_atom(x_+i*ca+j*cb, y_+i*sa+j*sb, z_, th_,fi_,at);
}
void cluster::delete_layer(numero x_, numero y_, numero z_,
                           numero a, numero b, numero al, numero be)
{
  int i,j,imax,jmax;
  numero d, ca=a*cos(al),cb=b*cos(al+be), sa=a*sin(al),sb=b*sin(al+be);
  if((Rmax/infinity)>0.9)  on_error(foutput,"delete_layer","Rmax too large");
  d=Rmax + sqrt(sqr(x_-x0)+sqr(y_-y0));
  imax=int(floor(d/a+1));  jmax=int(floor(d/b+1));
  for(i=-imax; i<=imax; i++)
  for(j=-jmax; j<=jmax; j++)  delete_atom(x_+i*ca+j*cb, y_+i*sa+j*sb, z_);
}
void cluster::add_surface(numero x_, numero y_, numero z_, numero a, numero c,
                          char *name, numero th_, numero fi_, int at)
{
  int i,imax, hcp=0;
  numero d, dx,dy,dz,ax,ay,bx,by, a1,a2, al,be, s2=sqrt(2), s3=sqrt(3);
  if((Rmax/infinity)>0.9)  on_error(foutput,"add_surface","Rmax too large");
  if(!strcmpC(name,"sc100")) {                                
    dx=0;          dy=0;         dz=-a;
    ax=a;          ay=0;         bx=0;          by=a;
  }  else
  if(!strcmpC(name,"sc110")) {                                
    dx=a/s2;       dy=0;         dz=-a/s2;
    ax=a*s2;       ay=0;         bx=0;          by=a;
  }  else
  if(!strcmpC(name,"sc111")) {                                
    dx=a*s2/s3;    dy=0;         dz=-a/s3;
    ax=a*s3/s2;    ay=-a/s2;     bx=a*s3/s2;    by=a/s2;
  }  else
  if(!strcmpC(name,"bcc100")) {                               
    dx=a/2;        dy=a/2;       dz=-a/2;
    ax=a;          ay=0;         bx=0;          by=a;
  }  else
  if(!strcmpC(name,"bcc110")) {                               
    dx=a/s2;       dy=0;         dz=-a/s2;
    ax=a/s2;       ay=-a/2;      bx=a/s2;       by=a/2;
  }  else
  if(!strcmpC(name,"bcc111")) {                               
    dx=a*s2/s3;    dy=0;         dz=-a/2/s3;
    ax=a*s3/s2;    ay=-a/s2;     bx=a*s3/s2;    by=a/s2;
  }  else
  if(!strcmpC(name,"fcc100")) {                               
    dx=a/2;        dy=0;         dz=-a/2;
    ax=a/2;        ay=-a/2;      bx=a/2;        by=a/2;
  }  else
  if(!strcmpC(name,"fcc110")) {                               
    dx=a/2/s2;     dy=a/2;       dz=-a/2/s2;
    ax=a/s2;       ay=0;         bx=0;          by=a;
  }  else
  if(!strcmpC(name,"fcc111")) {                               
    dx=a/s2/s3;    dy=0;         dz=-a/s3;
    ax=a*s3/2/s2;  ay=-a/2/s2;   bx=a*s3/2/s2;  by=a/2/s2;
  }  else
  if(!strcmpC(name,"hcp0001")) {                              
    dx=a/s3;       dy=0;         dz=-c;
    ax=a*s3/2;     ay=-a/2;      bx=a*s3/2;     by=a/2;
    hcp=1;
  }  else  on_error(foutput, "add_surface", "undefined surface type:", name);
  a1=sqrt(sqr(ax)+sqr(ay));
  a2=sqrt(sqr(bx)+sqr(by));
  al=acos(ax/a1);   if(ay<0) al=-al;                          
  be=ABS(acos((ax*bx+ay*by)/(a1*a2)));
  if(ax*bx+ay*by<0)  be=pi-be;
  d=Rmax + ABS(z_-z0);
  imax=int(floor(ABS(d/dz)+1));
    for(i=0; i<=imax; i++)
      if(hcp==0)
        add_layer(x_+i*dx, y_+i*dy, z_+i*dz, a1,a2, al,be, th_,fi_,at);
      else
        add_layer(x_+(i%2)*dx, y_+(i%2)*dy, z_+i*dz, a1,a2, al,be, th_,fi_,at);
}
void cluster::add_half_layer(numero x_, numero y_, numero z_,
                             numero a, numero b,
                             numero al, numero be, numero ff,
                             numero th_, numero fi_, int at)
{
  int i,j,imax,jmax;
  numero d, ca=a*cos(al),cb=b*cos(al+be), sa=a*sin(al),sb=b*sin(al+be);
  numero xx, yy, cf=cos(ff), sf=sin(ff);
  if((Rmax/infinity)>0.9)  on_error(foutput,"add_half_layer",
                                    "Rmax too large");
  d=Rmax + sqrt(sqr(x_-x0)+sqr(y_-y0));
  imax=int(floor(d/a+1));  jmax=int(floor(d/b+1));
  for(i=-imax; i<=imax; i++)
  for(j=-jmax; j<=jmax; j++) {
    xx=x_+i*ca+j*cb;
    yy=y_+i*sa+j*sb;
    if(cf*(xx-x_)+sf*(yy-y_)>-eps_d)  add_atom(xx, yy, z_, th_,fi_,at);
} }
void cluster::add_row(numero x_, numero y_, numero z_, numero a,
                      numero tt, numero ff, numero th_, numero fi_, int at)
{
  int i,imax;
  numero d, ax=a*sin(tt)*cos(ff), ay=a*sin(tt)*sin(ff), az=a*cos(tt);
  if((Rmax/infinity)>0.9)  on_error(foutput,"add_row","Rmax too large");
  d=Rmax + sqrt(sqr(x_-x0)+sqr(y_-y0)+sqr(z_-z0));
  imax=int(floor(d/a+1));
  for(i=-imax; i<=imax; i++)
    add_atom(x_+i*ax, y_+i*ay, z_+i*az, th_,fi_,at);
}
void cluster::delete_row(numero x_, numero y_, numero z_,
                         numero a, numero tt, numero ff)
{
  int i,imax;
  numero d, ax=a*sin(tt)*cos(ff), ay=a*sin(tt)*sin(ff), az=a*cos(tt);
  if((Rmax/infinity)>0.9)  on_error(foutput,"delete_row","Rmax too large");
  d=Rmax + sqrt(sqr(x_-x0)+sqr(y_-y0)+sqr(z_-z0));
  imax=int(floor(d/a+1));
  for(i=-imax; i<=imax; i++)  delete_atom(x_+i*ax, y_+i*ay, z_+i*az);
}
void cluster::add_emitter(int a)             
{
  add_emitter(a, atom[a]);
}
void cluster::add_emitter(int a, int sc)     
{
  int *emit_,*emit_s_,add,i;
  for(i=0, add=1; i<n_emit && add; i++)  if(emit[i]==a)  add=0;
  if(add) {
    emit_=new int [n_emit+1];    emit_s_=new int [n_emit+1];
    for(i=0; i<n_emit; i++)  {emit_[i]=emit[i];  emit_s_[i]=emit_s[i];}
    if(n_emit>0)  {delete [] emit;  delete [] emit_s;}
    emit_[n_emit]=a;    emit_s_[n_emit]=sc;
    emit=emit_;    emit_s=emit_s_;    n_emit++;
} }
void cluster::add_emitter_plane(numero x_, numero y_, numero z_,
                                numero th_, numero fi_)
{
  int i;  numero cth=cos(th_),sth=sin(th_),cfi=cos(fi_),sfi=sin(fi_);
  for(i=0; i<n; i++)                         
    if(ABS( ( (x[i]-x_)*cfi
             +(y[i]-y_)*sfi)*sth
           +  (z[i]-z_)*cth) < eps_d)  add_emitter(i,atom[i]);
}
void cluster::add_emitter_column(numero x_, numero y_)
{
  int i;
  for(i=0; i<n; i++)
    if(sqrt(sqr(x[i]-x_)+sqr(y[i]-y_))<eps_d)  add_emitter(i,atom[i]);
}
void cluster::all_emitters(int ta, int sc)
{
  int a;
  for(a=0; a<nn; a++)  if(atom[a]==ta)  add_emitter(a,sc);
}
int cluster::delete_non_emitters(void)
{
  numero *x_,*y_,*z_,*th_,*fi_;
  int i,j, *atom_,*atZ_;
  if(n_emit<=0)  return 1;
  x_=new numero [n_emit];  y_=new numero [n_emit];  z_=new numero [n_emit];
  atom_=new int [n_emit];  atZ_=new int [n_emit];
  th_=new numero [n_emit];  fi_=new numero [n_emit];
  for(i=j=0; i<nn && j<n_emit; i++) if(is_emitter(i)) {
    x_[j]=x[i];  y_[j]=y[i];  z_[j]=z[i];
    atom_[j]=atom[i];  atZ_[j]=atZ[i];  th_[j]=th[i];  fi_[j]=fi[i];  j++;
  }
  free_coordinates();
  x=x_;  y=y_;  z=z_;  atom=atom_;  atZ=atZ_;  th=th_;  fi=fi_;
  nn=n_emit;  set_n();
  free_emitters();
  return 0;
}
void cluster::displace(numero xx, numero yy, numero zz)
{
  int i;
  for(i=0; i<nn; i++) {x[i]+=xx;  y[i]+=yy;  z[i]+=zz;}
}
void cluster::rotate(double &xx, double &yy, double ang)
{
  numero rho,phi;
  rho=sqrt(xx*xx+yy*yy);
  if(!xx&&!yy) phi=0; else phi=atan2(yy,xx);
  xx=rho*cos(phi+ang);
  yy=rho*sin(phi+ang);
}
void cluster::rotate(int axis, numero ang)
{
  numero mux,muy,muz, rr;
  int i;
  for(i=0; i<nn; i++) {
    if(axis==0)  rotate(y[i],z[i],ang);  else    
    if(axis==1)  rotate(z[i],x[i],ang);  else    
    if(axis==2)  rotate(x[i],y[i],ang);          
    if(oriented) {                               
      mux=sin(th[i])*cos(fi[i]);
      muy=sin(th[i])*sin(fi[i]);
      muz=cos(th[i]);
      if(axis==0)  rotate(muy,muz,ang);  else    
      if(axis==1)  rotate(muz,mux,ang);  else    
      if(axis==2)  rotate(mux,muy,ang);          
      cartesian_to_spherical(mux,muy,muz, rr,th[i],fi[i]);
} } }
int cluster::is_emitter(int i)               
{                                            
  int j;                                     
                                             
  if(n_emit==0)  return 1;
  for(j=0; j<n_emit; j++)  if(emit[j]==i)  return 1;
  return 0;
}
void cluster::report_emitters(void)          
{                                            
  int i;                                     
  if(n_emit<1)  fprintf(foutput, "--- no emitters\n");               else
  if(n_emit==1) fprintf(foutput, "--- 1 emitter (in Angstroms)\n");  else
                fprintf(foutput, "--- %d emitters (in Angstroms)\n",n_emit);
  for(i=0; i<n_emit; i++)
    fprintf(foutput, "   %3d  %3d  %15.9f %15.9f  %15.9f\n",
            i+1, atom[emit[i]],
            x[emit[i]]*a0_au, y[emit[i]]*a0_au, z[emit[i]]*a0_au);
}
multiple_scattering::multiple_scattering(void)
{
  method=0;                 
  eta=1;                    
}
void multiple_scattering::free(void)
{                           
  o.free();  n.free();
}
void multiple_scattering::clear(void)
{                           
  o.clear();  n.clear();
}
void multiple_scattering::alloc(void)
{                           
  o.alloc(coord.n);
  n.alloc(coord.n);
}
void multiple_scattering::exact(void)
{
  int i,c,l,m, j,cp,lp,mp, n_c=o.g[0].n_c;
  numero *ri, *ii, *rj, *ij;
  complex *a, val;
  for(i=0; i<coord.n; i++) {
    n.g[i].free();
    n.g[i].alloc(scat.lmax_aa,n_c);
  }
  for(i=0, a=scat.aa.p;  i<coord.n;  i++)
  for(c=0, ri=n.g[i].re, ii=n.g[i].im;  c<n_c;  c++)
  for(l=0;  l<=scat.lmax_aa;  l++)
  for(m=-l;  m<=l;  m++, ri++, ii++)
  for(j=0;  j<coord.n;  j++)
  for(cp=0, rj=o.g[j].re, ij=o.g[j].im;  cp<n_c;  cp++)
  for(lp=0;  lp<=scat.lmax_aa;  lp++)
  for(mp=-lp;  mp<=lp;  mp++, a++)  if(lp<=o.g[j].lmax)  {
    val= (*a) * complex(*rj,*ij);
    *ri += real(val);
    *ii += imag(val);
    rj++;  ij++;
} }
void multiple_scattering::next_scattering(void)
{                                           
  int a,b;  vector gg;                      
  wave_function np;                         
                                            
  np.alloc(coord.n);
                                            
  for(a=0; a<coord.n; a++) {                
    np.g[a].alloc(-1,n.g[a].n_c);           
                                            
    scat.get_t_matrix(coord.atom[a]);
    scat.nn_c=n.g[a].n_c;
    if(scat.t_mat->ta>=0) {                 
                                            
      for(b=0; b<coord.n; b++)              
      if(scat.get_ij(a,b)) {                
        if((method==1 || method==2) && b<a)  gg.copy(np.g[b]); 
        else                                 gg.copy(n.g[b]);  
        if(gg.lmax>=0) {                    
          scat.get_ij();
          scat.propagate(gg, scat.t_mat->lmax);  
          np.g[a].add(gg);
        }
        gg.free();                          
        time_update();                      
      }
      scat.R_1TR(np.g[a],a,scat.t_mat);     
      if(method!=-1)  np.g[a].add(o.g[a]);  
    }  else  if(method!=-1)  np.g[a].add(o.g[a]);
  }
  for(a=0; a<coord.n; a++)
    if(method==2)  {n.g[a].ext(1-eta);  np.g[a].ext(eta);
                    n.g[a].add(np.g[a]);  np.g[a].free();}   
    else           {n.g[a].free();  n.g[a].copy(np.g[a]);}   
  np.free();
}
void multiple_scattering::write(FILE *fout, int o_, int max_order_)
{
  int  a,c,l,m;
  numero  *re,*im;
  if(o_==0) {
    fprintf(fout,"Momentum (a.u.)   %g  %g\n", real(scat.k), imag(scat.k));
    fprintf(fout,"Maximum order     %d\n", max_order_);
    fprintf(fout,"lmax_eff          %d\n\n", lmax_eff);
    fprintf(fout,"Cluster\n\n");            
    coord.write(fout,"l(A)");
    fprintf(fout,"\n");
    fprintf(fout,"Coefficients\n");
  }
  for(a=0; a<coord.n; a++) {
    fprintf(fout,"\n");
    fprintf(fout,"order %d   atom %d   nc %d   lmax %2d\n",
            o_, a+1, n.g[a].n_c, n.g[a].lmax);
    if(n.g[a].lmax>0)  fprintf(fout,"\n");
    re=n.g[a].re;  im=n.g[a].im;
    for(c=0; c<n.g[a].n_c; c++)
    for(l=0; l<=n.g[a].lmax; l++)
    for(m=-l; m<=l; m++, re++, im++)
      fprintf(fout,"  %1d %2d %3d   %16g %16g\n", c+1,l,m, *re,*im);
  }
  fprintf(fout,"\n");
}
int multiple_scattering::read(FILE *fin, int o_)
{
  char name[str_length_comments];
  int  a,c,l,m, n_c,lmax;
  numero  *re,*im;
  float  x,y;
  if(o_==0) {
    fscanf(fin,"%s%s%f%f", name, name, &x,&y);
    fscanf(fin,"%s%s%d", name, name, &max_order);
    fscanf(fin,"%s", name);   free_lmax();
    fscanf(fin,"%d", &lmax);
    if(lmax>lmax_max)
      on_error(foutput, "wf input", "lmax is too large in wf");
    fscanf(fin,"%s", name);   coord.read(fin);
    scat.set_lmax();                      
    if(lmax>lmax_eff)  init_lmax(lmax);
    final.init_Ylm();
    scat.set_k(complex(x,y));
    calc.kk[0]=sqr(real(scat.k))/2*au_eV;
    fscanf(fin,"%s", name);
    o.free();  o.alloc(coord.n);
    n.free();  n.alloc(coord.n);
  }  else  n.clear();
 
  for(a=0; a<coord.n; a++) {
    fscanf(fin,"%s%d", name, &o_);
    fscanf(fin,"%s%d", name, &a);    a--;
    fscanf(fin,"%s%d", name, &n_c);
    fscanf(fin,"%s%d", name, &lmax);
    
    n.g[a].alloc(lmax,n_c);
    re=n.g[a].re;  im=n.g[a].im;
    for(c=0; c<n.g[a].n_c; c++)
    for(l=0; l<=n.g[a].lmax; l++)
    for(m=-l; m<=l; m++, re++, im++) {
      fscanf(fin,"%d%d%d%f%f\n", &c,&l,&m, &x,&y);  c--;
      *re=x;  *im=y;
    }
    if(o_==0) o.g[a].copy(n.g[a]);
  }
  if(o_==max_order)  return 1;  else  return 0;
}
propagation::propagation(void)
{
  n_d=n_th=n_fi=n_t=n_G=nn_d=nn_th=nn_fi=nn_t=nn_G=0;  nn_c=1;
  dmax=infinity;
  translation_type=0;
  attenuation_type=0;
  t_mat0=NULL;
  iimfp.n=0;
  iimfp_flag=0;            
  T=0;
  int i; for(i=0; i<n_Debye_temperature_elements; i++) T_Debye[i]=100;
  V0=Vi=0;
  screening_length=1/a0_au;
  scattering_so=0;
  separation=1;
  Gab_cg_flag=0;
  n_a=0;
  a0=NULL;
  epsilon=0.00005;
  dim=0;   
  aa.n=aa.m=aa.nm=0;
  aa_output_flag=0;
}
void propagation::R_1TR(vector &gg, int i, t_matrix *tt)
{
  if(coord.oriented)  if(coord.th[i]!=0) {
    get_th(coord.th[i]);   ifi=get_fi(coord.fi[i]);   Rab(gg);  
    gg.spin_rotation(0,  coord.th[i], pi-coord.fi[i]);          
  }
  gg.tmat(tt);
  if(coord.oriented)  if(coord.th[i]!=0) {            Rab_1(gg);
    gg.spin_rotation(coord.fi[i]-pi, -coord.th[i], 0);
} }
void propagation::R_1TR(vector &gg, int i)
{
  get_t_matrix(coord.atom[i]);
  R_1TR(gg, i, t_mat);
}
void propagation::Rab(vector &gg)
{
  if(ifi<0) gg.Rfi(Rfi,2);         
  else      gg.Rfi(Rfi,0);         
            gg.Rth(Rth,0);         
}
void propagation::Rab_1(vector &gg)
{
            gg.Rth(Rth,1);         
  if(ifi<0) gg.Rfi(Rfi,3);         
  else      gg.Rfi(Rfi,1);         
}
void propagation::propagate(vector &gg, int lmax_final)
{
  if(separation>0) {
    Rab(gg);                     
    gg.Gz(Gz, lmax_final);       
    Rab_1(gg);                   
  }  else  Gab->prod(gg);
}
int propagation::get_ij(int i, int j)
{
  numero x,y,z;              
                             
  if(i==j)  x=y=z=d=th=fi=0; 
  else {
    x=coord.x[i]-coord.x[j]; 
    y=coord.y[i]-coord.y[j]; 
    z=coord.z[i]-coord.z[j]; 
    cartesian_to_spherical(x,y,z, d,th,fi);
    if(d>dmax || (coord.surface && coord.zs-coord.z[i]+d>dmax))  return 0;
  }
                                         
  if(i==j && separation>=0)  return 0;   
                                         
  return 1;
}
void propagation::get_ij(void)
{                       
  if(separation>0) {    
        get_d();        
        get_th(th);     
    ifi=get_fi(fi);     
  }  else               
  get_Gab();            
}
int propagation::get_d(void)
{                                 
  int n;                          
                                  
  if(n_d) {
    Gz=Gz0;  n=n_d;
    while(n--)
      if(fabs(d-Gz->d)<coord.eps_d                    
         && Gz->trans==translation_type)  return 0;   
      else  Gz=Gz->next;                              
  }  else  Gz0=NULL;                                  
  Gz=new green_z;  n_d++;  nn_d=n_d;  nn_G=0;         
  Gz->next=Gz0;   Gz0=Gz;                             
  Gz->init(k,d,translation_type,attenuation_type);
  return 0;
}
int propagation::get_th(numero th_)
{                                 
  int n, hz,l,m1,m2,v2l1;         
  rotation_th *Rtp;
  numero *pval,*pnval,*pth,*ppth,*th_save;
  if(n_th) {
    Rth=Rth0;  n=n_th;
    while(n--)  if(fabs(th_-Rth->th)<coord.eps_th)  return 0;
                else  Rth=Rth->next;
  }  else  Rth0=NULL;
  Rth=new rotation_th;  n_th++;  nn_th=n_th;
  Rth->next=Rth0;   Rth0=Rth;   Rth->th=th_;
  if(fabs(th_)   <coord.eps_th)  {Rth->val[0]=2;  return 0;}  
  if(fabs(th_-pi)<coord.eps_th)  {Rth->val[0]=3;  return 0;}  
  for(Rtp=Rth->next, hz=1, n=n_th; hz && --n;)
    if(fabs(pi-Rtp->th-th_)<coord.eps_th)  hz=0;  else Rtp=Rtp->next;
  if(hz) {Rth->init(th_); return 0;}
  th_save=new numero [sqr(2*lmax_eff+1)];
  for(l=0, pval=Rtp->val,pnval=Rth->val; l<=lmax_eff; l++) {
    v2l1=2*l+1;                                   
    for(m1=-l,pth=th_save;  m1<=l; m1++)          
    for(m2=-l; m2<=l; m2++, pval++,pth++)         
      *pth=*pval;                                 
    for(m1=-l,pth=th_save+v2l1;  m1<=l; m1++)     
    for(m2=-l,ppth=(--pth); m2<=l; m2++, pnval++,ppth+=v2l1)
      *pnval=((l&1)^(m1&1))?-(*ppth):(*ppth);
  }
  delete [] th_save;
  return 0;
}
int propagation::get_fi(numero fi_)
{                                 
  int n;                          
                                  
  if(n_fi) {                      
    Rfi=Rfi0;  n=n_fi;
    while(n--)  if(fabs(fi_-Rfi->fi)<coord.eps_th)    return  1;  else
                if(fabs(fi_-pi-Rfi->fi)<coord.eps_th ||
                   fabs(fi_+pi-Rfi->fi)<coord.eps_th) return -1;  else
                Rfi=Rfi->next;
  }  else  Rfi0=NULL;
  Rfi=new rotation_fi;  n_fi++;  nn_fi=n_fi;
  Rfi->next=Rfi0;   Rfi0=Rfi;   Rfi->fi=fi_;
  Rfi->init(pi-fi_);                          
  return 1;
}
int propagation::get_t_matrix(int ta)
{
  int n;                          
                                  
  if(n_t) {                       
    t_mat=t_mat0;  n=n_t;
    while(n--)  if((t_mat->ta)==ta)  return 0;  else  t_mat=t_mat->next;
  }  else  t_mat0=NULL;
  t_mat=new t_matrix;  n_t++;  nn_t=n_t;
  t_mat->next=t_mat0;  t_mat0=t_mat;  t_mat->ta=ta;
  if(ta>=0)  {get_scatterer(ta);  a->init_t_matrix(t_mat);}
  return 1;
}
void propagation::free_Green_functions(void)
{                                              
  while(n_d) {  n_d--;                         
    Gz=Gz0->next;   Gz0->free();
    delete Gz0;     Gz0=Gz;
  }
  if(Gab_cg_flag) {delete [] Gab_cg;  Gab_cg_flag=0;}
  while(n_G) {  n_G--;                         
    Gab=Gab0->next;   Gab0->free();            
    delete Gab0;      Gab0=Gab;
} }
void propagation::free_rotations(void)
{                                              
  while(n_th) {
    Rth=Rth0->next;   Rth0->free();   delete Rth0;
    Rth0=Rth;   n_th--;
  }
  while(n_fi) {
    Rfi=Rfi0->next;   Rfi0->free();   delete Rfi0;
    Rfi0=Rfi;   n_fi--;
} }
void propagation::free_t_matrix(void)
{                                              
  while(n_t) {
    t_mat=t_mat0->next;   t_mat0->free();   delete t_mat0;
    t_mat0=t_mat;   n_t--;
} }
int propagation::get_scatterer(int ta, int opt)
{
  int n=n_a;                      
                                  
  if(ta<0)  return 0;             
  a=a0;
  while(n--)  if((a->ta)==ta)  return 0;  else  a=a->next;
  if(opt) on_error(foutput,"get_scatterer","scatterer",ta,"undefined");
  return 1;
}
int propagation::free_scatterer(int ta)
{
  int n;  atomic_scattering *previous;
  if(n_a>0) {                     
    if(a0->ta==ta) {
      a=a0->next;  a0->free();  delete a0;  n_a--;  a0=a;
      return 0;
    }
    if(n_a>1) {
      a=a0->next;
      previous=a0;
      n=n_a-1;
      while(n--)
        if((a->ta)==ta) {
          n_a--;  previous->next=a->next;  a->free();  delete a;
          return 0;
        }  else  {previous=a;  a=a->next;}
  } }
  on_warning(foutput,"free_scatterer","scatterer",ta,"undefined");
  return 0;
}
void propagation::free_scatterers(void)
{
  while(n_a) {                               
    a=a0->next;  a0->free();  delete a0;     
    a0=a;  n_a--;
} }
void propagation::add_scatterer(int ta)
{
  int i;                                     
                                             
  a=a0;  i=1;
  while(a!=NULL && i)                        
    if(a->ta==ta) {                          
      free_scatterer(ta);                    
      i=0;                                   
    } else  a=a->next;                       
                                             
  a=new atomic_scattering;                   
  a->ta=ta;                                  
  a->next=a0;                                
  a0=a;   n_a++;                             
}                                            
void propagation::add_screened_scatterer(int ta, int ta_old)
{
  int n,Z,i, magnetic;
  numero *r,*V,*V2, rr;
  get_scatterer(ta_old);
  n=a->n_rad;  Z=a->atZ;  r=a->rad;  V=a->re;  V2=a->im;  magnetic=a->magnetic;
  add_scatterer(ta);
  a->select=3;  a->lmax=infinite_int;  a->n_rad=n;  a->atZ=Z;
  a->magnetic=magnetic;
  a->rad=copy(r,n);  a->re=copy(V,n);  if(magnetic) a->im=copy(V2,n);
  for(i=0; i<n; i++) {
    rr=a->rad[i];
    a->re[i] += -exp(-rr/screening_length)/rr;
    if(magnetic)  a->im[i] += -exp(-rr/screening_length)/rr;
} }
void propagation::scan_scatterer(int ta, char *name, int sel)
{
  FILE *fout=NULL;                           
                                             
  if(ta>=0) get_scatterer(ta);               
                                             
  if(strcmpC(name,"none"))  fout=open_file(foutput,name,"w");
  if(calc.verbose && sel!=2 && sel!=3) {
    fprintf(foutput,"--- properties of scatterer %d\n", ta);
    if(ta>=0) a->scan(foutput, sel);
  }
  if(fout!=NULL) {
    if(sel!=2 && sel!=3)  fprintf(fout,"--- properties of scatterer %d\n", ta);
    if(ta>=0) a->scan(fout, sel);
    close_file(fout);
} }
void propagation::set_k(int ik)
{
  numero kr,ki;
  if(calc.k_flag==1) {
    if(iimfp_flag==2)  set_k(sqrt(sqr(calc.k[ik])+2*V0+2*i_c*Vi));
    else {
      kr=sqrt(sqr(calc.k[ik])+2*V0);            
      if(iimfp_flag==0)  ki=iimfp.val(kr)/2;    
      else               ki=iimfp_TPP(kr)/2;    
      set_k(complex(kr,ki));
  } }  else  if(calc.k_flag==2)  set_k(calc.kc[ik]);
}
void propagation::set_k(complex k_)
{
  k=k_;
  refraction=V0/(sqr(real(k))/2-V0);            
  final.init_refraction(refraction);
}
numero propagation::iimfp_TPP(numero kr)
{
  numero E=sqr(kr)/2*au_eV;
  numero beta=-0.10+0.944/sqrt(TPP_Ep*TPP_Ep+TPP_Eg*TPP_Eg)
                   +0.069*pow(TPP_rho,0.1);
  numero gamma=0.191/sqrt(TPP_rho);
  numero U=TPP_Ep*TPP_Ep/829.4;
  numero C=1.97-0.91*U;
  numero D=53.4-20.8*U;
  numero imfp=E/(TPP_Ep*TPP_Ep*(beta*log(gamma*E)-C/E+D/(E*E)))/a0_au;
  return  1/imfp;
}
void propagation::scan_imfp(char *name)
{
  FILE *fout=NULL;                           
  if(strcmpC(name,"none"))  fout=open_file(foutput,name,"w");
  if(calc.verbose)  scan_imfp(foutput);
  if(fout!=NULL)   {scan_imfp(fout);  close_file(fout);}
}
void propagation::scan_imfp(FILE *fout)
{
  int ik;
  fprintf(fout,"--- electron inelastic mean free path\n");
  fprintf(fout,"     %5.5s  imfp(Angs)\n", str_units(calc.k_units));
  for (ik=0; ik<calc.n_k; ik++) {
    set_k(ik);
    fprintf(fout, " %9.4f", calc.kk[ik]);
    fprintf(fout, "   %9.4f\n", 1/(2*imag(k))*a0_au);
} }
void propagation::read_imfp(FILE *fin, char *my_file)
{
  FILE *fimfp;  int j, n_k;  numero k_,lam;
  char units1[str_length_units];                      
  char units2[str_length_units];    numero factor2;   
  if(!strcmpC(my_file,"inline"))  fimfp=fin;
  else                            fimfp=open_file(foutput,my_file,"r");
  n_k=read_int(fimfp);
  if(n_k>1)  read_name(fimfp,units1);
  read_name(fimfp,units2);
  factor2=units_conversion_factor(units2);
  iimfp.alloc(n_k,interpolation_type);
  if(n_k>1) {
    for(j=0; j<n_k; j++) {
      k_=read_numero(fimfp);
      k_=units_k(units1,k_);
      lam=read_numero(fimfp);
      if(lam<=0)  on_error(foutput,"read_imfp","negative i.m.f.p. =",lam,"");
      iimfp.put(j,k_,1/(factor2*lam));
    }
    iimfp.init();
  }  else {
    lam=read_numero(fimfp);
    if(lam<=0)  on_error(foutput,"read_imfp","negative i.m.f.p. =",lam,"");
    iimfp.init(1/(factor2*lam),0);
  }
  if(strcmpC(my_file,"inline"))  fclose(fimfp);
}
void propagation::set_muffin_tin_potentials(void)
{
  Mufpot s;
  int i,j, nn;
  free_scatterers();
  coord.sort(-1);                   
  for(i=0; i<coord.nn; i++)  coord.atom[i]=coord.atZ[i];
  s.init(coord.n, coord.nn, coord.atom, coord.x, coord.y ,coord.z);
  for(j=0; j<s.NR; j++) {
    add_scatterer(j);               
    nn=s.JRMT[j];                   
    a->select=3;                    
    a->lmax=infinite_int;           
    a->atZ=s.atZ[j];                
    element_properties(s.atZ[j]);
    a->mass=element_mass*amu_au;    
    a->n_rad=nn;                    
    a->rad=new numero [nn];
    a->re=new numero [nn];
    for(i=0; i<nn; i++) {
      a->rad[i]=s.RX[i];            
      a->re[i]=s.V[j][i];           
  } }
  s.free();
}
void propagation::set_lmax(void)
{
  int lmax=-1, i;
  for(i=0; i<coord.nn; i++)
    if(coord.atom[i]>=0) {                
      get_scatterer(coord.atom[i]);       
      if((a->lmax)>lmax)  lmax=a->lmax;   
    }                                     
  init_lmax(lmax);
  if(lmax_eff<0)  on_error(foutput,"set_lmax","lmax =",lmax_eff,"< 0");
}
int propagation::get_n_c(void)
{
  int i;
  for(i=0; i<coord.nn; i++)  if(coord.atom[i]>=0) {
      get_scatterer(coord.atom[i]);
      if((a->select==0 && a->n_c==2)     ||   
         (a->select==1 && a->n_c==2)     ||   
         (a->select==2 && a->n_c==2)     ||   
         (a->select==3 && scattering_so) ||   
         (a->select==3 && a->magnetic)   ||   
          a->select==4)  return 2;
    }
  return 1;                                   
}
void propagation::report_links(void)
{
  int n_at=calc.n_atoms_last;
  if(nn_G>0)
  fprintf(foutput,"              Number of Gab's           = %3d\n", nn_G);
  else {
  fprintf(foutput,"--- number of elements used in the last calculation\n");
  fprintf(foutput,"              Number of atoms           = %3d\n", n_at);
  fprintf(foutput,"              Number of bond distances  = %3d\n", nn_d);
  fprintf(foutput,"              Number of polar angles    = %3d\n", nn_th);
  fprintf(foutput,"              Number of azim. angles    = %3d\n", nn_fi);
  fprintf(foutput,"              Number of t matrices      = %3d\n", nn_t);
  fprintf(foutput,"              Number of spin components = %3d\n", nn_c);
} }
void propagation::get_exact(int n_c, int lmax_)
{                                       
  int i,j, c,l,m, cp,lp,mp, ia,ja;      
  vector gg;                            
  numero *rr,*ii;                       
                                        
  if(lmax_>lmax_eff)  lmax_=lmax_eff;   
  lmax_aa=lmax_;                        
  dim=coord.n*n_c*lp12[lmax_aa];
  aa.alloc(dim,dim);
  for(i=0; i<coord.n; i++) {
    get_t_matrix(coord.atom[i]);
    for(j=0; j<coord.n; j++)  if(get_ij(i,j)) {
      for(cp=0; cp<n_c; cp++)
      for(lp=0; lp<=lmax_aa; lp++)
      for(mp=-lp; mp<=lp; mp++) {
        gg.free();
        gg.alloc(lp,n_c);
        gg.re[cp*lp12[lp]+lp*lp+lp+mp]=1;
        get_ij();
        if(t_mat->lmax<lmax_aa)  propagate(gg, t_mat->lmax);
        else                     propagate(gg, lmax_aa);
        R_1TR(gg, i, t_mat);
        for(c=0, rr=gg.re, ii=gg.im; c<n_c; c++)
        for(l=0; l<=lmax_aa; l++)
        for(m=-l; m<=l; m++)  if(l<=gg.lmax)  {
          ia=(i*n_c + c) *lp12[lmax_aa] +l*l   +l  +m;
          ja=(j*n_c + cp)*lp12[lmax_aa] +lp*lp +lp +mp;
          aa.a(ia,ja, -complex(*rr,*ii));      
          rr++;  ii++;
  } } } }
  for(i=0; i<dim; i++)  aa.add(i,i,1);         
  --aa;   gg.free();                           
  FILE *fout_aa;
  if(aa_output_flag) {
    fout_aa=open_file(foutput, aa_output, "a");
    fprintf(fout_aa, "%d %d %d %g\n",
            coord.n, n_c, lmax_aa, real(sqr(k))/2*au_eV);
    for(i=0; i<coord.n; i++)
      fprintf(fout_aa, "%d %g %g %g\n", i, coord.x[i],coord.y[i],coord.z[i]);
    for(i=0; i<dim; i++)
    for(j=0; j<dim; j++)
      fprintf(fout_aa, "%g %g\n", real(aa.a(i,j)), imag(aa.a(i,j)));
    close_file(fout_aa);
} }
void propagation::free_exact(void)              
{                                              
  if(dim>0) {aa.free();  dim=0;}
}
int propagation::get_Gab(void)
{                                 
  int n;                          
                                  
  if(n_G) {
    Gab=Gab0;  n=n_G;
    while(n--)
      if(fabs(d -Gab->d) <coord.eps_d   &&            
         fabs(th-Gab->th)<coord.eps_th  &&            
         fabs(fi-Gab->fi)<coord.eps_th)   return 0;   
      else  Gab=Gab->next;                            
  }  else  Gab0=NULL;
  Gab=new green;  n_G++;  nn_G=n_G;                   
  Gab->next=Gab0;   Gab0=Gab;                         
  Gab->init(k,d,th,fi, separation);                   
                                                      
  return 0;                                           
}
atomic_scattering::atomic_scattering(void)
{
  n_k=0;                 
  comp=0;                
  select=0;              
  n_rad=0;  rad=NULL;    
  lmax=-1;  n_rad=0;  atZ=0;
  mass=infinite;
  magnetic=0;
}
void atomic_scattering::free(void)
{
  int l;
  if(select==1 || select==2) {
    if(n_k>0)  {delete [] re;  delete [] im;  delete [] k;  n_k=0;}
  }  else
  if(select==3) {
    if(n_rad>0) {delete [] rad;  delete [] re;  if(magnetic) delete [] im;}
    n_rad=0;  rad=NULL; atZ=0;  magnetic=0;
  }  else
  if(lmax>=0 && select!=-1) {
    for(l=0; l<=lmax; l++)  {rr[l].free();  if(comp) ii[l].free();}
    delete [] rr;
    if(comp)  delete [] ii;
    lmax=-1;
} }
int atomic_scattering::read(FILE *fin, char *name)
{
  char units[str_length_units],             
       idx[str_length_comments];            
  FILE *ffin;                               
  free();
  if(!strcmpC(name,"inline"))  ffin=fin;
  else                         ffin=open_file(foutput,name,"r");
  read_name(ffin,idx);
  if(!strcmpC(idx,"potential") || !strcmpC(idx,"potential2")) {
    if(n_rad>0) {delete [] rad;  delete [] re; if(magnetic) delete [] im;}
    if(!strcmpC(idx,"potential2"))  magnetic=1;  else  magnetic=0;
    select=3;  lmax=infinite_int;  read_potential(ffin);  return 0;
  }
  n_k=read_int(idx);                        
  if(n_k>1)  read_name(ffin, units);        
  lmax=read_int(ffin);                      
                                            
  read_name(ffin, idx);                     
  if(!strcmpC(idx,"regular")){select=0;n_c=1; read_phases(ffin,units);} else
  if(!strcmpC(idx,"ex"))     {select=0;n_c=2; read_phases(ffin,units);} else
  if(!strcmpC(idx,"so"))     {select=4;n_c=2; read_phases(ffin,units);} else
  if(!strcmpC(idx,"axial"))  {select=1;       read_full  (ffin,units);} else
  if(!strcmpC(idx,"full"))   {select=2;       read_full  (ffin,units);} else
  on_error(foutput,"scatterer","wrong argument:",idx);
  if(strcmpC(name,"inline"))  fclose(ffin);
  return 0;
}
void atomic_scattering::read_potential(FILE *fin)
{
  char l_units[str_length_units], E_units[str_length_units];
  char name[str_length_comments];
  int  i, rV;
  n_rad=read_int(fin);  atZ=read_int(fin);
  read_name(fin,l_units);
  read_name(fin,E_units);
  read_name(fin,name);  if(!strcmpC(name,"rV"))  rV=1;  else  rV=0;
  rad=new numero [n_rad];
  re=new numero [n_rad];       
  if(magnetic)
  im=new numero [n_rad];       
                               
  for(i=0; i<n_rad; i++) {
    rad[i]=units_k(l_units, read_numero(fin));
    if(rV) re[i]=units_k(l_units, units_E(E_units, read_numero(fin)))/rad[i];
    else   re[i]=units_E(E_units, read_numero(fin));
    if(magnetic)
    if(rV) im[i]=units_k(l_units, units_E(E_units, read_numero(fin)))/rad[i];
    else   im[i]=units_E(E_units, read_numero(fin));
} }
void atomic_scattering::read_phases(FILE *fin, char *units)
{
  int j,l,p,c;      numero k_, x,y;
  char cmpl[str_length_comments];            
  read_name(fin, cmpl);                            
  if(!strcmpC(cmpl,"real"))     comp=0;  else      
  if(!strcmpC(cmpl,"complex"))  comp=1;  else      
                                comp=2;
            rr = new spline [n_c*(lmax+1)];
  if(comp)  ii = new spline [n_c*(lmax+1)];
  for(c=p=0; c<n_c; c++)
  for(l=0; l<=lmax; l++, p++)  {
              rr[p].alloc(n_k,interpolation_type);
    if(comp)  ii[p].alloc(n_k,interpolation_type);
  }
  for(j=0; j<n_k; j++) {                            
    if(n_k>1)  k_=units_k(units,read_numero(fin));  
    for(c=p=0; c<n_c; c++)
    for(l=0; l<=lmax; l++, p++) {
      x=read_numero(fin);
      if(comp)  y=read_numero(fin);  else  y=0;
      if(n_k==1)  {rr[p].init(x,0);   if(comp) ii[p].init(y,0);}
      else        {rr[p].put(j,k_,x); if(comp) ii[p].put(j,k_,y);}
      if((comp==0 || comp==1) && j>0) { 
        x=rr[p].y[j-1];  y=rr[p].y[j];
        while(rr[p].y[j]>x+pi)  rr[p].y[j]-=2*pi;
        while(rr[p].y[j]<=x-pi) rr[p].y[j]+=2*pi;
  } } }
  if(n_k>1)
  for(c=p=0; c<n_c; c++)
  for(l=0; l<=lmax; l++, p++) {              
    if(comp) ii[p].init();  rr[p].init();
} }
void atomic_scattering::read_full(FILE *fin, char *units)
{
  int i,j,n;          
  numero *r_,*i_;     
                      
  n_c=read_int(fin);
  if(select==1)  n=(sqr(lmax+1)+(lmax*(lmax+1)*(2*lmax+1))/3);  
  else           n=sqr(sqr(lmax+1));                            
  n=n*sqr(n_c);                                                 
  r_=re=new numero [n_k*n];
  i_=im=new numero [n_k*n];  k=new numero [n_k];
  for(i=0; i<n_k; i++) {
    if(n_k>1)  k[i]=units_k(units, read_numero(fin));
    else       k[i]=0;
    for(j=0; j<n; j++) {
      if(n_c>1) read_int(fin);                      
                read_int(fin);  read_int(fin);      
      if(n_c>1) read_int(fin);                      
                read_int(fin);  read_int(fin);      
      *r_=read_numero(fin);   r_++;
      *i_=read_numero(fin);   i_++;
} } }
void atomic_scattering::init_t_matrix(t_matrix *t)
{
  t->free();
  t->lmax=(lmax>lmax_eff)?lmax_eff:lmax;   
  if(select== 0)  init_t_matrix_diag(t);               else
  if(select== 1)  init_t_matrix_full(real(scat.k),t);  else
  if(select== 2)  init_t_matrix_full(real(scat.k),t);  else
  if(select== 3)  init_t_matrix_potential(t);          else
  if(select== 4)  init_t_matrix_so(t);                 else
  if(select==-1)  init_t_matrix_sphere(t);             else
  on_error(foutput,"init_t_matrix","select =", select, "not implemented");
}
complex *atomic_scattering::tl_phases(int n_c_, int &lmax_)
{
  complex *tl;                 
  int     c,l,p;               
  numero  psre,psim, expdi;
  if(lmax_>lmax)  lmax_=lmax;
  tl=new complex [n_c_*(lmax_+1)];
  for(c=p=0; c<n_c_; c++)
  for(l=0; l<=lmax_; l++, p++) {
              psre=rr[c*(lmax+1)+l].val(real(scat.k));
    if(comp)  psim=ii[c*(lmax+1)+l].val(real(scat.k));
    if(comp==0)  expdi=1;
    if(comp==1)  expdi=exp(-2*psim);
    if(comp==0 || comp==1)
          tl[p]=complex(sin(2*psre)*expdi/2, (1-cos(2*psre)*expdi)/2);
    else  tl[p]=complex(psre, psim);
  }
  cut_tl_phases(tl, n_c_, lmax_);
  if(select==4)  tl=tl_temperature_so(tl, lmax_);
  else           tl=tl_temperature(tl, n_c_, lmax_);
  return tl;
}
void atomic_scattering::cut_tl_phases(complex *tl, int n_c_, int &lmax_)
{
  int c,l,p, llmax=-1;     
  for(c=p=0; c<n_c_; c++)
  for(l=0; l<=lmax_; l++, p++) {
    if(scat.epsilon<mod(tl[p]) && l>llmax) llmax=l;
  }
  if(llmax<lmax_) {
    if(n_c_>1)  for(l=0; l<=llmax; l++)  tl[l+(llmax+1)]=tl[l+(lmax_+1)];
    lmax_=llmax;
} }
complex *atomic_scattering::tl_temperature(complex *tl, int n_c_, int &lmax_)
{
  if(scat.T<1e-10)  return tl;
  complex *ttl;  ttl=new complex [n_c_*(lmax_max+1)];
  numero  alpha=3*scat.T/(2*mass*kb_au*sqr(scat.T_Debye[atZ]));
  complex  alphak=2*alpha*mod2(scat.k), jl;
  int c,l,p, ll,lpp;
  for(c=p=0; c<n_c_; c++)
  for(l=0; l<=lmax_max; l++, p++)  ttl[p]=0;
  for(ll=0; ll<=lmax_max+lmax_; ll++) {
    jl=i_l(ll) * exp(-alphak) * besselj(ll,-i_c*alphak);
    for(c=p=0; c<n_c_; c++)
    for(l=0; l<=lmax_max; l++, p++)
    for(lpp=ABS(l-ll); lpp<=ABS(l+ll) && lpp<=lmax_; lpp+=2)
        ttl[p] += jl * tl[lpp+c*(lmax_+1)]
                  * sqrt(4*pi*(2*ll+1)*(2*lpp+1)/(2*l+1.0))
                  * Gaunt(lpp,ll,l,0,0,0);
  }
  lmax_=lmax_max;  cut_tl_phases(ttl, n_c_, lmax_);
  delete [] tl;  return ttl;
}
numero atomic_scattering::Gaunt_j(int l, int j, int mu, int lp, int jp,
                                  int lpp)
{
  if(ABS(j-l-l)>1 || ABS(jp-lp-lp)>1 || ABS(mu)>j)  return 0;
  numero val0=0, val;
  int s,mup,mpp;
  for(mup=-jp; mup<=jp; mup+=2) {
    mpp=(mu-mup)/2;
    if(ABS(mpp)<=lpp) {
      val=0;
      for(s=-1; s<=1; s+=2)
        val += CGspin1_2(l+l,j,s,mu) * CGspin1_2(lp+lp,jp,s,mup)
               * Gaunt(lp,lpp,l,(mup-s)/2,mpp,(mu-s)/2);
      val0+=sqr(val);
  } }
  return val0;
}
complex *atomic_scattering::tl_temperature_so(complex *tl, int &lmax_)
{
  if(scat.T<1e-10)  return tl;
  complex *ttl;  ttl=new complex [2*(lmax_max+1)];
  numero  alpha=3*scat.T/(2*mass*kb_au*sqr(scat.T_Debye[atZ]));
  complex  alphak=2*alpha*mod2(scat.k), jl;
  int l,j,c,p, lp,jp,cp, mu,lpp;
  for(c=p=0; c<2; c++)
  for(l=0; l<=lmax_max; l++, p++)  ttl[p]=0;
  for(lpp=0; lpp<=lmax_max+lmax_; lpp++) {
    jl=4*pi*i_l(lpp) * exp(-alphak) * besselj(lpp,-i_c*alphak);
    for(c=p=0; c<2; c++)
    for(l=0; l<=lmax_max; l++, p++) if(l>0 || c==0) {
      if(c==0) j=l+l+1; else j=l+l-1;
      for(lp=0; lp<=lmax_; lp++)
      for(cp=0; cp<2; cp++) if(lp>0 || cp==0) {
        if(cp==0) jp=lp+lp+1; else jp=lp+lp-1;
        ttl[p] += jl * tl[lp+cp*(lmax_+1)] * Gaunt_j(l,j,j,lp,jp,lpp);
  } } }
  lmax_=lmax_max;  cut_tl_phases(ttl, 2, lmax_);
  delete [] tl;  return ttl;
}
void atomic_scattering::init_t_matrix_diag(complex *tl, t_matrix *t)
{
  int c,l,p;                                  
                                              
  t->select=0;                                
  t->re = new numero [(t->n_c)*(t->lmax+1)];
  t->im = new numero [(t->n_c)*(t->lmax+1)];
  for(c=p=0; c<t->n_c; c++)
  for(l=0; l<=t->lmax; l++, p++) {
    t->re[p]=real(tl[p]);  t->im[p]=imag(tl[p]);
} }
void atomic_scattering::init_t_matrix_diag(t_matrix *t)
{
  complex *tl;
  t->n_c=n_c;  tl=tl_phases(t->n_c, t->lmax);
  init_t_matrix_diag(tl, t);
  delete [] tl;
}
void atomic_scattering::init_t_matrix_so(complex *tl, t_matrix *t)
{
  numero trp,tip,trm,tim;                     
  int l,m,p,q;
  t->select=4;
  t->re = new numero [4*lp12[t->lmax]];
  t->im = new numero [4*lp12[t->lmax]];
  for(l=p=q=0; l<=t->lmax; l++)
  for(m=-l; m<=l; m++, p+=8, q+=4)  {
    trp=real(tl[l]);  tip=imag(tl[l]);                
    if(l==0) trm=tim=0;                               
    else    {trm=real(tl[l+(t->lmax+1)]);  tim=imag(tl[l+(t->lmax+1)]);}
    t->re[q  ]= cgcg12[p  ]*trp + cgcg12[p+1]*trm;    
    t->re[q+1]= cgcg12[p+4]*trp + cgcg12[p+5]*trm;    
    t->re[q+2]= cgcg12[p+2]*trp + cgcg12[p+3]*trm;    
    t->re[q+3]= cgcg12[p+6]*trp + cgcg12[p+7]*trm;    
    t->im[q  ]= cgcg12[p  ]*tip + cgcg12[p+1]*tim;    
    t->im[q+1]= cgcg12[p+4]*tip + cgcg12[p+5]*tim;    
    t->im[q+2]= cgcg12[p+2]*tip + cgcg12[p+3]*tim;    
    t->im[q+3]= cgcg12[p+6]*tip + cgcg12[p+7]*tim;    
} }
void atomic_scattering::init_t_matrix_so(t_matrix *t)
{
  complex *tl;
  t->n_c=2;  tl=tl_phases(t->n_c, t->lmax);   
  init_t_matrix_so(tl, t);
  delete [] tl;
}
complex *atomic_scattering::tl_potential(int n_c_, int &lmax_)
{
  complex *tl;
  numero E=sqr(real(scat.k))/2, *ph, cc;
  Salvat Dirac;                             
  int l, ll, rel;
  if(magnetic)  ll=lmax_+1;  else  ll=n_c_*(lmax_+1);
  tl=new complex [n_c_*(lmax_+1)];
  if(n_rad<1)  on_error(foutput,"scatterer", "undefined MT potential");
  if(n_c_==1 || magnetic)  {cc=c_au; c_au*=1000; rel=0;}  else  rel=1;
  ph=new numero [ll];
    Dirac.phases(rad,re,n_rad, atZ, E, rel, lmax_, ph);
    for(l=0; l<ll; l++)  tl[l]=sin(ph[l])*exp(i_c*ph[l]);
  if(magnetic) {
    Dirac.phases(rad,im,n_rad, atZ, E, rel, lmax_, ph);
    for(l=0; l<ll; l++)  tl[l+ll]=sin(ph[l])*exp(i_c*ph[l]);
  }
  if(n_c_==1)  c_au=cc;
  delete [] ph;  Dirac.free();
  cut_tl_phases(tl, n_c_, lmax_);
  if(n_c_==1 || magnetic)  tl=tl_temperature(tl, n_c_, lmax_);
  else                     tl=tl_temperature_so(tl, lmax_);
  return tl;
}
void atomic_scattering::init_t_matrix_potential(t_matrix *t)
{
  complex *tl;
  if(scat.scattering_so || magnetic)  t->n_c=2;
  else                                t->n_c=1;
  tl=tl_potential(t->n_c, t->lmax);    
  if(scat.scattering_so)  init_t_matrix_so(tl, t);
  else                    init_t_matrix_diag(tl, t);
  delete [] tl;
}
complex *atomic_scattering::tl_sphere(int lmax_)
{
  complex *tl;
  complex r=scat.k*a;
  complex rp=sqrt(sqr(scat.k)+2*V)*a;
  complex jl,jlp;
  int     l;
  tl=new complex [lmax_+1];
  for(l=0; l<=lmax_; l++) {
    jl=besselj(l,rp);
    jlp=rp*besseljp(l,rp);
    tl[l]= (-jlp * besselj(l,r) + jl * r*besseljp(l,r))
           /(jlp * besselh(l,r) - jl * r*besselhp(l,r));
  }
  cut_tl_phases(tl, 1, lmax_);
  tl=tl_temperature(tl, 1, lmax_);
  return tl;
}
void atomic_scattering::init_t_matrix_sphere(t_matrix *t)
{
  complex *tl;
  t->n_c=1;  tl=tl_sphere(t->lmax);           
  init_t_matrix_diag(tl, t);
  delete [] tl;
}
void atomic_scattering::init_t_matrix_full(numero kk, t_matrix *t)
{
  int i,n, c,l,m,cp,lp,mp;
  numero x, *re_,*im_, *r1,*i1,*r2,*i2;
  if(select==1)  n=(sqr(lmax+1)+(lmax*(lmax+1)*(2*lmax+1))/3);  
  else           n=sqr(sqr(lmax+1));                            
  n=n*sqr(n_c);                                                 
  t->select=select;
  t->n_c=n_c;
  if(n_k==1) {t->re=copy(re,n);  t->im=copy(im,n);}
  else {
    if(kk<k[0] || k[n_k-1]<kk)
      on_error(foutput,"init_t_matrix_full", "momentum", kk, "out of range");
    i=1;  while(kk>k[i]) i++;
    x=(kk-k[i-1])/(k[i]-k[i-1]);
    t->re=copy_(re+(i-1)*n, re+i*n, n, x);
    t->im=copy_(im+(i-1)*n, im+i*n, n, x);
  }
  if(t->lmax<lmax) {                                
    r1=t->re;  i1=t->im;
    if(select==1)  n=sqr(t->lmax+1)+(t->lmax*(t->lmax+1)*(2*t->lmax+1))/3;
    else           n=sqr(sqr(t->lmax+1));      n=n*sqr(n_c);
    r2=re_=new numero [n];
    i2=im_=new numero [n];
    for(c=0; c<n_c; c++)
    for(l=0; l<=lmax; l++)
    for(m=-l; m<=l; m++)
    for(cp=0; cp<n_c; cp++)
    for(lp=0; lp<=lmax; lp++)
    for(mp=-lp; mp<=lp; mp++)  if(select==2 || m==mp) {
      if(l<=t->lmax && lp<=t->lmax) {*r2=*r1;  *i2=*i1;  r2++;  i2++;}
      r1++;  i1++;
    }
    delete [] t->re;  delete [] t->im;
    t->re=re_;  t->im=im_;
} }
numero *atomic_scattering::copy_(numero *a_, numero *b_, int n, numero x)
{
  int i;
  numero *c;
  c=new numero [n];
  for(i=0; i<n; i++)  c[i]=(1-x)*a_[i]+x*b_[i];
  return c;
}
int atomic_scattering::scan(FILE *fout, int sel)
{
  numero k_,th, dp,dm, mu,Pl0,Pl1, ps;
  complex f,g, val, *tl;
  int ik,j,l, so_on, ex_on, n_c_, lmax_;
  if(sel==-1) {
    if(select!=-1 && select!=3 && n_k<1)
    fprintf(fout, "    no scattering properties are defined\n"); else {
    if(select==-1)  fprintf(fout, "    sphere of constant potential\n");
    if(select== 0)  fprintf(fout, "    non-relativistic phase shifts\n");
    if(select== 1)  fprintf(fout, "    t matrix for axially-symmetric atom\n");
    if(select== 2)  fprintf(fout, "    t matrix for arbitrary atom\n");
    if(select== 3)  fprintf(fout, "    muffin-tin potential\n");
    if(select== 4)  fprintf(fout, "    relativistic phase shifts\n");}
    fprintf(fout, "    atomic mass = ");
    if(mass<infinite/2)  fprintf(fout, "%g a.u.\n", mass);
    else                 fprintf(fout, "infinite\n");
    if(select==3) fprintf(fout, "    atomic number = %d\n", atZ);
    if(select!=-1 && select!=3) {
      fprintf(fout, "    lmax = %d\n", lmax);
      fprintf(fout, "    number of momenta = %d\n", n_k);
      if(n_k>1)  fprintf(fout, "    momentum in the range %g-%g a.u.\n",
                               k[0], k[n_k-1]);
    }
    return 0;
  }
  if(sel==3) {
    if(select!=3)  on_warning(fout,"scan scatterer potential",
                              "the muffin-tin potential is not stored");
    if(magnetic)  fprintf(fout, "potential2\n");
    else          fprintf(fout, "potential\n");
    fprintf(fout, "%d %d l(au) E(au) rV\n", n_rad, atZ);
    for(j=0; j<n_rad; j++)  {
       fprintf(fout, "%g %g", rad[j], rad[j]*re[j]);
       if(magnetic)  fprintf(fout, " %g\n", rad[j]*im[j]);
       else          fprintf(fout, "\n");
    }
    return 0;
  }
  if(select==1 || select==2) {
    on_warning(fout, "scan scatterer", "scatterer defined by its t matrix");
    return 0;
  }
  if(select==0 && n_k<1) {
    on_warning(fout, "scan scatterer", "empty set of phase shifts");
    return 0;
  }
  init_lmax(lmax_max);
  if(select==-1 || select==0 || (select==3 && scat.scattering_so==0))
        so_on=0;
  else  so_on=1;
  if((select==0 && n_c==2) || (select==3 && magnetic))  ex_on=1;
  else                                                  ex_on=0;
  if(scat.scattering_so || (select==3 && magnetic))  n_c_=2;  else  n_c_=1;
  if(sel==0) {
    if(calc.n_k>1)    fprintf(fout, "     eV");
    if(final.n_th>1)  fprintf(fout, "  theta");
    fprintf(fout, "   cross sect.");
    if(so_on==0)  fprintf(fout, "        Re{f}        Im{f}\n");
    else          fprintf(fout, "            S      S |f|^2\n");
    for(ik=0; ik<calc.n_k; ik++) {
      scat.k=k_=calc.k[ik];  lmax_=lmax_max;
      if(select==-1)  tl=tl_sphere(lmax_);           else
      if(select== 3)  tl=tl_potential(n_c_, lmax_);  else
      if(select== 0)  tl=tl_phases(n_c, lmax_);      else
      if(select== 4)  tl=tl_phases(2, lmax_);        else
      on_error(foutput, "scan scatterer", "undefined option");
      for(j=0; j<final.n_th; j++) {
        th=final.th_out[j];  mu=cos(th);
        f=g=0;
        
        for(l=0; l<=lmax_; l++) {
          if(l && so_on) Pl1=legendre(l,1,mu);
          Pl0=legendre(l,0,mu);
          if(so_on) {dp=tl_deltal(tl[l]);  dm=tl_deltal(tl[l+lmax_+1]);}
	  else
          if(select==0 && n_c==2) dp=tl_deltal((tl[l]+tl[l+lmax_+1])/2.0);
          else                    dp=tl_deltal(tl[l]);
          if(so_on) {                   
            f+=1/(2*i_c*k_)
               *((l+1)*(exp(2*i_c*dm)-1.0)+l*(exp(2*i_c*dp)-1.0))*Pl0;
            if(l)  g+=1/(2*i_c*k_)*(exp(2*i_c*dp)-exp(2*i_c*dm))*Pl1;
          }  else
            f+=1/k_*(2*l+1)*tl[l]*Pl0;  
        }
        
        if(calc.n_k>1)    fprintf(fout, " %6.2f", au_eV*sqr(k_)/2);
        if(final.n_th>1)  fprintf(fout, " %6.2f", 180/pi*th);
        if(so_on)
          fprintf(fout, "  %12g %12g %12g\n",
                  mod2(f)+mod2(g),                        
                  -2*imag(f*conj(g))/(mod2(f) + mod2(g)), 
                  -2*imag(f*conj(g))                      
                 );
        else
          fprintf(fout, "  %12g %12g %12g\n",
                  mod2(f),                                
                  real(f),                                
                  imag(f)                                 
                 );
      }
      delete [] tl;
  } }
  if(sel==2) {
    fprintf(fout, "%d", calc.n_k);
    if(calc.n_k>1)  fprintf(fout, " %s", calc.k_units);
    fprintf(fout, " %d", lmax_max);
    if(so_on)  fprintf(fout, " so");  else
    if(ex_on)  fprintf(fout, " ex");  else
               fprintf(fout, " regular");
    fprintf(fout, " tl\n");
    for(ik=0; ik<calc.n_k; ik++) {
      scat.k=k_=calc.k[ik];  lmax_=lmax_max;
      if(select==-1)  tl=tl_sphere(lmax_);           else
      if(select== 3)  tl=tl_potential(n_c_, lmax_);  else
      if(select== 0)  tl=tl_phases(n_c, lmax_);      else
      if(select== 4)  tl=tl_phases(2, lmax_);        else
      on_error(foutput, "scan scatterer", "undefined option");
      if(calc.n_k>1)  fprintf(fout, "%6.2f", calc.kk[ik]);
      for(l=0; l<=lmax_max; l++)
        if(l<=lmax_)  fprintf(fout, " %g %g", real(tl[l]), imag(tl[l]));
        else   fprintf(fout, " 0 0");
      if(so_on || ex_on || (select==3 && magnetic))
      for(l=0; l<=lmax_max; l++)
        if(l<=lmax_)
          fprintf(fout, " %g %g", real(tl[l+lmax_+1]), imag(tl[l+lmax_+1]));
        else   fprintf(fout, " 0 0");
      fprintf(fout, "\n");
      delete [] tl;
  } }
  free_lmax();   return 0;
}
final_state::final_state(void)
{
  n_th=n_fi=0;
  n_1=n_2=0;
  Ylm0_th_flag=Ylm0_fi_flag=0;
  mesh_flag=0;
}
void final_state::init_th(numero thi, numero thf, int nth)
{
  int i;
  free_mesh();
  if(n_th>0) {delete [] th;  delete [] th_out;  delete [] transmission;}
  if(nth<1)  n_th=0;  else {
    n_th=nth;
    th=new numero [n_th];                 
    th_out=new numero [n_th];             
    transmission=new numero [n_th];       
    if(n_th==1)  {th[0]=th_out[0]=thi; transmission[0]=1;}  else
    for(i=0; i<n_th; i++) {
      th[i]=th_out[i]=thi+i*(thf-thi)/(n_th-1);
      transmission[i]=1;
} } }
void final_state::init_phi(numero fii, numero fif, int nfi)
{
  int j;
  free_mesh();
  if(n_fi>0) {delete [] fi;}
  if(nfi<1)  n_fi=0;  else {
    n_fi=nfi;
    fi=new numero [n_fi];                   
    if(n_fi==1)  fi[0]=fii;
    else  for(j=0; j<n_fi; j++)  fi[j]=fii+j*(fif-fii)/(n_fi-1);
} }
void final_state::init_refraction(numero refraction)
{
  int i;
  for(i=0; i<n_th; i++)
    if(refraction==0)   th[i]=th_out[i];
    else  {
      th[i]=asin(sin(th_out[i])/sqrt(1+refraction));
      if(th_out[i]>pi/2)  th[i]=pi-th[i];
    }
  init_transmission(refraction);
}
void final_state::init_transmission(numero refraction)
{
  int i;
  for(i=0; i<n_th; i++)
    if(refraction!=0)
          transmission[i]=cos(th_out[i])/cos(th[i]) / (1+refraction);
    else  transmission[i]=1;
}
void final_state::init_Ylm(void)
{                                           
  int i,j,m;                                
  if(n_th*n_fi) {
    free_Ylm();
    Ylm0_th_flag=n_th;                      
    Ylm0_fi_flag=n_fi;                      
    Ylm0_th=new array [n_th];               
    Ylm0_fi=new array [n_fi];               
                                            
    for(i=0; i<n_th; i++) {                 
      Ylm0_th[i].alloc(lp12[lmax_eff]);     
      init_Ylm(th[i],0,Ylm0_th[i].a);       
    }
    for(j=0; j<n_fi; j++) {                 
      Ylm0_fi[j].alloc(2*lmax_eff+1);
      for(m=-lmax_eff; m<=lmax_eff; m++)
        Ylm0_fi[j].a[m+lmax_eff]=exp(i_c*(m*fi[j]));
} } }
void final_state::init_Ylm(numero th_, numero fi_, complex *Ylm_)
{
  int l,m,lm;                               
  numero cth=cos(th_);                      
  for(l=lm=0; l<=lmax_eff; l++)
  for(m=-l; m<=l; m++,lm++)  Ylm_[lm]=Ylmmu(l,m,cth,fi_);
}
complex *final_state::Ylm(int i, int j)
{                                           
  int l,m,lm;                               
  complex *Ylm_;                            
  Ylm_=new complex [lp12[lmax_eff]];
  for(l=lm=0; l<=lmax_eff; l++)
  for(m=-l; m<=l; m++,lm++)  Ylm_[lm]=Ylm0_th[i].a[lm] *
                                      Ylm0_fi[j].a[m+lmax_eff];
  return Ylm_;
}
vector final_state::Ylm(complex *Ylm_)
{                                           
  numero *tr,*ti;                           
  vector temp;                              
  int l,m,lm;                               
  temp.alloc(lmax_eff,1);
  for(l=lm=0, tr=temp.re, ti=temp.im; l<=lmax_eff; l++)
  for(m=-l; m<=l; m++, lm++, tr++, ti++) {
    *tr=real(Ylm_[lm]);
    *ti=imag(Ylm_[lm]);
  }
  return temp;
}
void final_state::init_mesh(int n_1_, int n_2_)
{
  free();  n_1=n_1_;  n_2=n_2_;
  x=new numero [n_1*n_2];
  y=new numero [n_1*n_2];
  z=new numero [n_1*n_2];
}
void final_state::init_mesh(numero x0, numero y0, numero z0,
                            numero x1, numero y1, numero z1,
                            numero x2, numero y2, numero z2,
                            int n_1_, int n_2_)
{
  init_mesh(n_1_,n_2_);
  numero h1,h2;  int i,j;
  if(n_1>1)  h1=1/(n_1-1.0);  else  h1=0;
  if(n_2>1)  h2=1/(n_2-1.0);  else  h2=0;
  mesh_flag=1;               
  for(i=0; i<n_1; i++)
  for(j=0; j<n_2; j++) {
    x[i*n_2+j]=x0+i*(x1-x0)*h1+j*(x2-x1)*h2;
    y[i*n_2+j]=y0+i*(y1-y0)*h1+j*(y2-y1)*h2;
    z[i*n_2+j]=z0+i*(z1-z0)*h1+j*(z2-z1)*h2;
} }
void final_state::init_mesh(numero r,
                            numero th1, numero th2,
                            numero fi1, numero fi2, int n_1_, int n_2_)
{
  init_mesh(n_1_,n_2_);
  numero tth,ffi, cth,sth, h1,h2;
  if(n_1>1)  h1=1/(n_1-1.0);  else  h1=0;
  if(n_2>1)  h2=1/(n_2-1.0);  else  h2=0;
  int i,j;
  mesh_flag=2;               
  for(i=0; i<n_1; i++) {
    tth=th1+i*(th2-th1)*h1;
    cth=cos(tth);  sth=sin(tth);
    for(j=0; j<n_2; j++) {
      ffi=fi1+i*(fi2-fi1)*h2;
      x[i*n_2+j]=r*sth*cos(ffi)+coord.x0;
      y[i*n_2+j]=r*sth*sin(ffi)+coord.y0;
      z[i*n_2+j]=r*cth+coord.z0;
} } }
void final_state::free(void)
{
  if(n_th>0) {
    n_th=0;  delete [] th;  delete [] th_out;  delete [] transmission;
  }
  if(n_fi>0) {n_fi=0;  delete [] fi;}
  free_mesh();
  free_Ylm();
}
void final_state::free_Ylm(void)            
{                                           
  int i;                                    
                                            
  if(Ylm0_th_flag)  {
    for(i=0; i<Ylm0_th_flag; i++)  Ylm0_th[i].free();
    delete [] Ylm0_th;  Ylm0_th_flag=0;
  }
  if(Ylm0_fi_flag)  {
    for(i=0; i<Ylm0_fi_flag; i++)  Ylm0_fi[i].free();
    delete [] Ylm0_fi;  Ylm0_fi_flag=0;
} }
void final_state::free_mesh(void)
{
  if(mesh_flag>0) {
    delete [] x;  delete [] y;  delete [] z;
    mesh_flag=n_1=n_2=0;
} }
#ifdef jga_erf            
#else                     
#define jga_erf 1         
                          
                          
numero gammln(numero xx)
{
  numero stp=2.5066282746310005,
         x=xx-1,tmp,ser=1.000000000190015,
         cof[6];
  int j;
  cof[0]=76.18009172947146;      cof[1]=-86.50532032941677;
  cof[2]=24.01409824083091;      cof[3]=-1.231739572450155;
  cof[4]=0.1208650973866179e-2;  cof[5]=-0.5395239384953e-5;
  x=xx;
  tmp=x+5.5;
  tmp=(x+0.5)*log(tmp)-tmp;
  for(j=0; j<6; j++)  {x++;  ser+=cof[j]/x;}
  return tmp+log(stp*ser/xx);
}
numero gamma(numero x)
{
  if(x<0)  return gamma(x+1)/x;
           return exp(gammln(x));
}
numero gser(numero a, numero x, int p)
{
  int n, itmax=101;
  numero sum,del,ap, eps=3.0e-7, gln=gammln(a);
  if(x<=0) {
    if(x<0) printf("Error in gser: x<0\n");
    return 0;
  } else {
    ap=a;    sum=1/a;    del=sum;
    for(n=1; n<itmax; n++) {
      ap++;
      del=del*x/ap;
      sum+=del;
      if(ABS(del)<ABS(sum)*eps) n=itmax+1;
    }
    if(n<=itmax)  printf("Error in gser: a too large, itmax too small\n");
    if(p)  return  sum*exp(-x+a*log(x)-gln);
    else   return  exp(gln)-sum*exp(-x+a*log(x));
} }
numero gcf(numero a, numero x, int p)
{
  int n, itmax=101;
  numero gold,g,fac,b1,b0,anf,ana,an,a1,a0, eps=3.0e-7, gln=gammln(a);
  gold=0;  a0=1;  a1=x;  b0=0;  b1=1;  fac=1;
  for(n=1; n<itmax; n++) {
    an=1*n;    ana=an-a;    a0=(a1+a0*ana)*fac;
    b0=(b1+b0*ana)*fac;     anf=an*fac;
    a1=x*a0+anf*a1;         b1=x*b0+anf*b1;
    if(a1!=0) {
      fac=1/a1;   g=b1*fac;
      if(ABS((g-gold)/g)<eps) n=itmax+1;
      gold=g;
  } }
  if(n<=itmax)  printf("Error in gcf: a too large, itmax too small\n");
  if(p)  return  exp(-x+a*log(x)-gln)*g;
  else   return  exp(-x+a*log(x))*g;
}
numero gammp(numero a, numero x)
{
  if(x<0 || a<=0) {printf("Error in gammp: invalid arguments\n");  return 0;}
  if(x<a+1)  return  gser(a,x,1);
  else       return  1-gcf(a,x,1);
}
numero gammq(numero a, numero x)
{
  if(x<0 || a<=0) {printf("Error in gammq: invalid arguments\n");  return 0;}
  if(x<a+1)  return  1-gser(a,x,1);
  else       return gcf(a,x,1);
}
numero gamma(numero a, numero x)
{
  if(a<0)  return (gamma(a+1,x)-pow(x,a)*exp(-x))/a;
  if(x==0)   return gamma(a);
  if(x<0) {printf("Error in gammq: invalid arguments\n");  return 0;}
  if(x<a+1)  return gser(a,x,0);
  else       return gcf(a,x,0);
}
numero erf(numero x)
{
  if(x<0) return -gammp(0.5,sqr(x));
  else    return  gammp(0.5,sqr(x));
}
numero erfc(numero x)
{
  if(x<0) return 1+gammp(0.5,sqr(x));
  else    return gammq(0.5,sqr(x));
}
#ifdef jga_complex
complex erf(complex xx)
{
  if(imag(xx)<0)   return -erf(-xx);  else
  if(imag(xx)==0)  return complex(erf(real(xx)),0);
  numero x=real(xx), y=imag(xx), ex,xy,x2,_2x, cxy,sxy, ey, ch,sh, fn,gn;
  complex val,z;  val=0;
  int n,n2,nmax=10;
  x2=x*x;  ex=exp(-x2);  xy=2*x*y;  cxy=cos(xy);  sxy=sin(xy);  _2x=2*x;
  for(n=1; n<=nmax; n++) {
     n2=n*n;  ey=exp(n*y);  ch=(ey+1/ey)/2;  sh=(ey-1/ey)/2;
     fn = _2x*(1-ch*cxy) + n*sh*sxy;
     gn = _2x*   ch*sxy  + n*sh*cxy;
     val=val+(exp(-n2/4.0)/(n2+4*x2))*complex(fn,gn);
  }
  if(x!=0)  z=complex(1-cxy,sxy)/_2x;  else  z=complex(0,y);
  return  erf(x) + (ex/pi) * (z + 2*val);
}
complex erfc(complex x) {return 1.0-erf(x);}
complex erfw_asymptotic(complex z)
{
  complex z2=z*z;
  if(real(z)>6 || imag(z)>6)
        return i_c*z*(0.5124242/(z2-0.2752551)+0.05176536/(z2-2.724745));
  else  return i_c*z*(0.4613135/(z2-0.1901635)+0.09999216/(z2-1.7844927)
                      + 0.002883894/(z2-5.5253437));
}
complex erfw_series(complex z, int n)
{
  complex fct1, fct2, ff=-z*z, val;
  int i;
  fct1=1;
  fct2=i_c*z/(sqrt(pi)/2);  val=fct1+fct2;
  for(i=1; i<=n; i++) {
    fct1=fct1*ff/(i+0.0);
    fct2=fct2*ff/(i+0.5);
    val+=fct1+fct2;
  }
  return val;
}
complex erfw(complex z)            
{
  numero x=real(z), y=imag(z);
  if(y<0)           return 2*exp(-z*z) - erfw(-z);
  if(x<0)           return conj(erfw(complex(-x,y)));
  if(x>3.9 || y>3)  return erfw_asymptotic(z);
                    return erfw_series(z, 100);
}
complex erfw(complex z, complex b)
{
  numero x=real(z), y=imag(z);
  complex b2=b*b;
  complex eb2=exp(-b2);
  if(y<0)           return 2*exp(-z*z-b2) - eb2 * erfw(-z);
  if(x<0)           return eb2 * conj(erfw(complex(-x,y)));
  if(x>3.9 || y>3)  return eb2 * erfw_asymptotic(z);
                    return eb2 * erfw_series(z, 100);
}
complex gamma1_2(int n, complex x)       
{
  if(n==0) return sqrt(pi)*exp(-x)*erfw(euler(sqrt(mod(x)),(arg(x)+pi)/2));
  if(n<0)  return (gamma1_2(n+1,x)-pow(x,n+0.5)*exp(-x))/(n+0.5);
           return (n-0.5)*gamma1_2(n-1,x)+pow(x,n-0.5)*exp(-x);
}
#endif  
#endif  
#ifdef jga_Kambe          
#else                     
#define jga_Kambe 1       
                          
                          
class Kambe {
  numero  a0,b0;                          
  numero  ax,ay, bx,by;                   
  numero  Gax,Gay, Gbx,Gby;               
  numero  th0,phi0;                       
  numero  qx,qy,qz;                       
  complex k;                              
  numero  att;                            
  numero  eps;                            
  numero  a,a_2;                          
  int D3lm(void);
  #define integl_ij 40401                 
  complex  integl(int ll, complex kRj);
  void     integl_clear(void);
  int     *integl_flag[integl_ij], integl_l, integl_i;
  complex *integl_data[integl_ij];
  numero   integl_x0, integl_y0;
  complex Almz(int l_, int m_,
               numero x, numero y, numero z, numero z1_, int flag_);
  void tested_sum(complex coef, int ij, int index);
  void perimeter(complex coef,
                 int ij0, int ij1, int index);   
  
  numero  x0,y0,z0, z1, area;
  int     l,m, i,j,s, n,lm, flag;
  complex k2, sum, qzc;
public:
  int ijmin, ijstep;                      
  Kambe(void) {ijmin=5; ijstep=3; integl_l=-1; x0=y0=0;}
  void free(void);
  void integl_alloc(int ll);
  void init(complex k_, numero th0_, numero phi0_, numero eps_);
  void init(complex k_, numero eps_) {init(k_,pi,0,eps_);}
  void lattice(numero a0_, numero b0_, numero al, numero be);
  void lattice(numero a0_) {lattice(a0_,a0_,0,pi/2);}
  void attenuation(numero att_) {att=att_;}
  void q(numero qx_, numero qy_, numero qz_) {qx=qx_;  qy=qy_;  qz=qz_;}
  complex Alm(int l, int m, numero x, numero y, numero z);
  complex Blm(int l, int m, numero x, numero y, numero z, numero z1_);
  complex Blm_NaCl1(int l, int m, int n_);
  complex Blm_NaCl2(int l, int m, int n_);
  complex Clm(int l, int m, numero x, numero y, numero z, numero z1_);
  complex Alm(int l, int m, numero x, numero y, numero z,
                            numero qx_, numero qy_, numero qz_);
  complex Almj(int l, int m, numero x, numero y, numero z, int j_);
  void report(void)   
  {
     printf("q=(%g,%g,%g)  att=%g\n", qx,qy,qz, att);
     printf("k=(%g,%g)\n", real(k), imag(k));
     printf("a=(%g,%g)  b=(%g,%g)\n", ax,ay, bx,by);
} };
void Kambe::lattice(numero a0_, numero b0_, numero al, numero be)
{
  a0=a0_;         b0=b0_;            
  ax=a0*cos(al);  bx=b0*cos(al+be);  
  ay=a0*sin(al);  by=b0*sin(al+be);  
  area=ax*by-ay*bx;                  
  Gax= 2*pi * by/area;               
  Gay=-2*pi * bx/area;               
  Gbx=-2*pi * ay/area;               
  Gby= 2*pi * ax/area;
  area=ABS(area);
}
void Kambe::init(complex k_, numero th0_, numero phi0_, numero eps_)
{
  init_ffact();  free();
  k=k_;  k2=k*k;
  th0=th0_;  phi0=phi0_;
  eps=eps_;
                                     
  a=mod2(k)*area/(32*pi);            
  if(a>8)  a=8;                      
  a_2=sqrt(a);                       
  att=imag(k)/ABS(cos(th0));
  qx=real(k)*sin(th0)*cos(phi0);
  qy=real(k)*sin(th0)*sin(phi0);
  qz=real(k)*cos(th0);
}
void Kambe::free(void)
{
  int ij;
  if(integl_l>=0)
    for(ij=0; ij<integl_ij; ij++) {
      delete [] integl_flag[ij];  delete [] integl_data[ij];
    }
  integl_l=-1;
}
complex Kambe::Alm(int l_, int m_, numero x, numero y, numero z)
{
  complex coef_D1lm, coef_D2lm;
  l=l_;  m=m_;  x0=x;  y0=y;  z0=z;
  if(mod(k*z0)>1e-10)  return Almz(l,m,x0,y0,z0,0,0);       
  if(!l && !m && x0==0 && y0==0)  sum=i_c*k/sqrt(4*pi);     
  else                            sum=0;
  if((l+m)%2==0) {    
    lm=(l-ABS(m))/2;
                                                            
    coef_D1lm=-1/(k*area) * i_l(1-m)/pow(2.0,l+0.0)         
              * sqrt((2*l+1)*fact[l+m]*fact[l-m]);
    coef_D2lm=-k/(4*pi)*(sign_1l(l+(l+m)/2)/pow(2.0,l+0.0)) 
              * sqrt((2*l+1)*fact[l-m]*fact[l+m])
              / (fact[(l-m)/2]*fact[(l+m)/2]);
    D3lm();                                                 
    perimeter(coef_D1lm, 0, ijmin, 1);                      
    if(mod(k*(integl_x0-x0))+mod(k*(integl_x0-x0))>1e-9)  integl_clear();
    if(x0!=0 || y0!=0)  perimeter(coef_D2lm, 0,ijmin, 2);   
    else                perimeter(coef_D2lm, 1,ijmin, 2);
    tested_sum(coef_D1lm, ijmin, 1);                        
    tested_sum(coef_D2lm, ijmin, 2);                        
  }
  return  -4*pi * exp(i_c*(qx*x0+qy*y0)) * sum;
}
complex Kambe::Blm(int l_, int m_, numero x, numero y, numero z, numero z1_)
{
  complex val=0;     
  z1_=ABS(z1_);
  if(mod(k*z1_)<=1e-10)  on_error("Kambe::Blm","z1=0");
  while(z>=0 || mod(k*z)<=1e-10)  {val+=Alm(l_,m_,x,y,z);  z-=z1_;}
  return val + Almz(l_,m_,x,y,z,-z1_,1);
}
complex Kambe::Blm_NaCl1(int l_, int m_, int n_)
{
  return   (Blm(l_,m_, a0/2,    0,           (n_/2)*a0,  a0)
          + Blm(l_,m_,    0, a0/2,           (n_/2)*a0,  a0)
          + Blm(l_,m_,    0,    0, -a0/2+((n_+1)/2)*a0,  a0)
          + Blm(l_,m_, a0/2, a0/2, -a0/2+((n_+1)/2)*a0,  a0))
        * exp(-n_*(a0/2)*att);
}
complex Kambe::Blm_NaCl2(int l_, int m_, int n_)
{
  return   (Blm(l_,m_, a0/2,    0, -a0/2+((n_+1)/2)*a0,  a0)
          + Blm(l_,m_,    0, a0/2, -a0/2+((n_+1)/2)*a0,  a0)
          + Blm(l_,m_,    0,    0,           (n_/2)*a0,  a0)
          + Blm(l_,m_, a0/2, a0/2,           (n_/2)*a0,  a0))
        * exp(-n_*(a0/2)*att);
}
complex Kambe::Clm(int l_, int m_, numero x, numero y, numero z, numero z1_)
{
  numero att_=att;   
  complex val;
  z1_=ABS(z1_);
  if(mod(k*z1_)<=1e-10)  on_error("Kambe::Clm","z1=0");
  while(z>=z1_)  z-=z1_;
  while(z<=-z1_) z+=z1_;
  attenuation(0);
  val=  Almz(l_,m_,x,y,z+z1_,z1_,1)
      + Almz(l_,m_,x,y,z-z1_,-z1_,1)
      + Alm (l_,m_,x,y,z);
  attenuation(att_);
  return val;
}
complex Kambe::Alm(int l_, int m_, numero x, numero y, numero z,
                   numero qx_, numero qy_, numero qz_)
{
  numero qqx=qx, qqy=qy, qqz=qz;   qx=qx_;  qy=qy_;  qz=qz_;
  complex  val=Alm(l,m,x,y,z);     qx=qqx;  qy=qqy;  qz=qqz;
  return val;
}
complex Kambe::Almj(int l_, int m_, numero x, numero y, numero z, int j_)
{
  numero  h=real(k)/10000;
  if(j_==0)  return -i_c*( Alm(l,m,x,y,z,qx+h,qy,qz)
                          -Alm(l,m,x,y,z,qx-h,qy,qz))/(2*h);  else
  if(j_==1)  return -i_c*( Alm(l,m,x,y,z,qx,qy+h,qz)
                          -Alm(l,m,x,y,z,qx,qy-h,qz))/(2*h);  else
  if(j_==2)  return -i_c*( Alm(l,m,x,y,z,qx,qy,qz+h)
                          -Alm(l,m,x,y,z,qx,qy,qz-h))/(2*h);
  return 0;
}
void Kambe::integl_clear(void)
{
  int ij,ll;
  if(integl_l>=0)  for(ij=0; ij<integl_ij; ij++)
                   for(ll=-1; ll<=integl_l; ll++)  integl_flag[ij][ll+1]=0;
  integl_x0=x0;  integl_y0=y0;  integl_i=0;
}
void Kambe::integl_alloc(int ll)
{
  int ij;  integl_l=ll;
  for(ij=0; ij<integl_ij; ij++) {
    integl_flag[ij]=new int [integl_l+2];
    integl_data[ij]=new complex [integl_l+2];
  }
  integl_clear();
}
complex Kambe::integl(int ll, complex kRj)                  
{
  if(integl_i<integl_ij && ll<=integl_l)
    if(integl_flag[integl_i][ll+1])   return integl_data[integl_i][ll+1];
  complex kR2a=sqr(kRj)/(4*a), val;
  if(ll==-1)  val= sqrt(pi)        * exp(a) / (2*i_c)
                     * (  erfw( a_2+i_c*kRj/(2*a_2), kRj/(2*a_2))
                        - erfw(-a_2+i_c*kRj/(2*a_2), kRj/(2*a_2)));  else
  if(ll== 0)  val= sqrt(pi)/(kRj/2.0) * exp(a) /  2.0
                     * (  erfw( a_2+i_c*kRj/(2*a_2), kRj/(2*a_2))
                        + erfw(-a_2+i_c*kRj/(2*a_2), kRj/(2*a_2)));  else
  if(ll>0)    val= ((ll-0.5)*integl(ll-1,kRj)
                             - integl(ll-2,kRj)
                             + exp(a-kR2a)/pow(a,ll-0.5))/(a*kR2a);  else
              val=  (ll+1.5)*integl(ll+1,kRj)
                             - integl(ll+2,kRj) * a*kR2a
                             + exp(a-kR2a)/pow(a,ll+1.5);
  if(integl_i<integl_ij && ll<=integl_l) {
    integl_flag[integl_i][ll+1]=1;  integl_data[integl_i][ll+1]=val;
  }
  return val;
}
int Kambe::D3lm(void)                                     
{
  if(l || m || x0!=0 || y0!=0)  return 0;
  numero ea=exp(a);
  sum += k/(2*pi) * (ea/a_2 + i_c*sqrt(pi)*(erfw(a_2)*ea-1.0));
  return 0;
}
complex Kambe::Almz(int l_, int m_,
                    numero x, numero y, numero z, numero z1_, int flag_)
{
  x0=x;  y0=y;  z0=z;  l=l_;  m=m_;            
  z1=z1_;  flag=flag_;                         
                                               
  if(z0*z1<=0 && flag)  on_error("Kambe::Almz","z0*z1<=0");
  s=(z0>0)?1:-1;
  qzc=qz-i_c*att;                              
  sum=0;
  perimeter(1, 0, 4*ijmin, 3);
  tested_sum(1, 4*ijmin, 3);
  return (4*pi) * 2*pi/area * sign_1l(l+m) * i_l(l) * sum;
}
void Kambe::tested_sum(complex coef, int ij, int index)
{
  complex sum0;  ij++;
  do {
    sum0=sum;
    perimeter(coef, ij,ij+ijstep-1,index);  ij+=ijstep;
  } while(mod((sum-sum0)/sum)>eps);
}
void Kambe::perimeter(complex coef, int ij0, int ij1, int index)
{
  complex  partial_sum=0, sum1, fct1,fct2, Gam,cth, qGz,qGz2, kR;
  numero   Rx,Ry, rRx,rRy,rR, Gx,Gy, qGx,qGy, qG2, phi;
  for(i=-ij1; i<=ij1; i++)                     
  for(j=-ij1; j<=ij1; j++) {
    if(ABS(i)<ij0)  if(j==-ij0+1)  j=ij0;
    if(index==1) {                             
      Gx=i*Gax+j*Gbx;  Gy=i*Gay+j*Gby;
      qGx=qx+Gx;       qGy=qy+Gy;
      qG2=sqr(qGx)+sqr(qGy);
      if(qG2!=0 || m==0) {
        qGz2=k2-qG2;  qGz=sqrt(qGz2);
        fct1=pow(sqrt(qG2)/k,ABS(m))*pow(qGz/k,lm+lm-1.0);
        fct2=qG2/qGz2;
        for(n=lm,sum1=0; n>=0; n--) {
          sum1 += fct1 / (fact[n]*fact[(l-m)/2-n]*fact[(l+m)/2-n])
                       * gamma1_2(-n,-a*qGz2/k2-i_c*1e-15);
          fct1=fct1*fct2;
        }
        if(qGx==0 && qGy==0)  phi=0;  else  phi=atan2(qGy,qGx);
        partial_sum += exp(-i_c*(qGx*x0+qGy*y0+m*phi)) * sum1;
    } }  else
    if(index==2) {                             
      Rx=i*ax+j*bx;       Ry=i*ay+j*by;
      rRx=x0+Rx;          rRy=y0+Ry;
      rR=sqrt(sqr(rRx)+sqr(rRy));   kR=k*rR;
      sum1=pow(kR/2.0,l) * integl(l,kR);      integl_i++;
      if(rRx==0 && rRy==0)  phi=0;  else  phi=atan2(-rRy,-rRx);
      partial_sum += exp(-i_c*(qx*Rx+qy*Ry+m*phi)) * sum1;
    }  else
    if(index==3) {                             
      Gx=i*Gax+j*Gbx;  Gy=i*Gay+j*Gby;
      qGx=-qx+Gx;      qGy=-qy+Gy;
      qG2=sqr(qGx)+sqr(qGy);
      Gam=sqrt(qG2-k2-complex(0,1e-15));
      if(qGx==0 && qGy==0)  phi=0;  else  phi=atan2(qGy,qGx);
      cth=s*i_c*Gam/k;
      sum1=  exp(-Gam*ABS(z0) + i_c*(Gx*x0+Gy*y0+qzc*z0))/Gam
           * Ylmmu(l,-m,cth,phi);
      if(flag)  sum1=sum1/(1.0-exp((-s*Gam+i_c*qzc)*z1));
      partial_sum+=sum1;
  } }
  sum +=coef*partial_sum;
}
#endif  
void green::init(complex k, numero dd, numero tth, numero ffi, int periodic)
{
  int c,l,m,cp,lp,mp,lpp,mpp, i;
  numero *Gr,*Gi, cth, jl,yl, *cg, krd,kid;
  complex *hL;   hL=new complex [sqr(2*lmax_eff+1)];
  complex hh, kd, val;
  Kambe a;
  if(scat.Gab_cg_flag==0) {
    for(l=i=0; l<=lmax_eff; l++)
    for(m=-l; m<=l; m++)
    for(lp=0; lp<=lmax_eff; lp++)
    for(mp=-lp; mp<=lp; mp++)
    for(lpp=ABS(l-lp); lpp<=l+lp; lpp+=2) i++;
    scat.Gab_cg=new numero [i];
    cg=scat.Gab_cg;
    for(l=0; l<=lmax_eff; l++)
    for(m=-l; m<=l; m++)
    for(lp=0; lp<=lmax_eff; lp++)
    for(mp=-lp; mp<=lp; mp++)
    for(lpp=ABS(l-lp); lpp<=l+lp; lpp+=2) {
      mpp=m-mp;
      *cg=Gaunt(lp,lpp,l,mp,mpp);
      cg++;
    }
    scat.Gab_cg_flag=1;
  }
  d=dd;  th=tth;  fi=ffi;  cth=cos(tth);
  kd=k*d;  krd=real(kd);  kid=imag(kd);
  re=new numero [sqr(lp12[lmax_eff])];  Gr=re;
  im=new numero [sqr(lp12[lmax_eff])];  Gi=im;
  if(periodic==0)         
    for(l=0; l<=2*lmax_eff; l++) {   
      if(scat.attenuation_type==0)  hh=besselilh(l,kd);
      else {
        besseljy(l,krd,jl,yl);
        hh=i_l(l) * (i_c*jl-yl) * exp(-kid);
      }
      hh=4*pi*hh;
      for(m=-l; m<=l; m++)  hL[l*l+l+m]= hh*conj(Ylmmu(l,m,cth,fi));
    }
  else {                  
    a.lattice(3.61/a0_au);
    a.init(k, 30*pi/180,50*pi/180, 0.0001);
    a.attenuation(0);
    for(l=0; l<=2*lmax_eff; l++)
    for(m=-l; m<=l; m++)
      hL[l*l+l+m]=a.Alm(l,m,d*sin(th)*cos(fi),d*sin(th)*sin(fi),d*cth) / k;
  }
  cg=scat.Gab_cg;
  for(c=0; c<1; c++)                           
  for(l=0; l<=lmax_eff; l++)
  for(m=-l; m<=l; m++)
  for(cp=0; cp<1; cp++)
  for(lp=0; lp<=lmax_eff; lp++)
  for(mp=-lp; mp<=lp; mp++, Gr++, Gi++) {
    val=0;
    mpp=m-mp;
    for(lpp=ABS(l-lp); lpp<=l+lp; lpp+=2) {
      val += hL[lpp*lpp+lpp+mpp] * (*cg);
      cg++;
    }
    *Gr=real(val);
    *Gi=imag(val);
  }
  delete [] hL;
}
void green::prod(vector &a)
{
  int c,l,m,cp,lp,mp;
  numero *Gr,*Gi,*ar,*ai,*br,*bi;
  vector b;
  b.copy(a);  a.free();  a.alloc(lmax_eff,b.n_c);
  Gr=re;  Gi=im;  ar=a.re;  ai=a.im;
  for(c=0; c<a.n_c; c++)                        
  for(l=0; l<=lmax_eff; l++)
  for(m=-l; m<=l; m++, ar++, ai++)
  for(cp=0, br=b.re, bi=b.im; cp<a.n_c; cp++)
  for(lp=0; lp<=lmax_eff; lp++)
  for(mp=-lp; mp<=lp; mp++, Gr++, Gi++)  if(lp<=b.lmax) {
      *ar += (*Gr)*(*br) - (*Gi)*(*bi);
      *ai += (*Gr)*(*bi) + (*Gi)*(*br);
      br++;  bi++;
  }
  b.free();
}
int electron_bin::init_wf_LEED(int spin_)
{
  int  a,l,m,lm, dim=lp12[lmax_eff];           
  complex ekRa,psi, *Ylm;                      
  vector gg;                                   
  numero kr=real(scat.k);
  numero fi0=calc.e_ffi;
  numero th0=calc.e_tth;
                                               
  k0x=kr*sin(th0)*cos(fi0);                    
  k0y=kr*sin(th0)*sin(fi0);                    
  k0z=kr*cos(th0);                             
                                               
  Ylm=new complex [dim];                       
  final.init_Ylm(th0,fi0,Ylm);                 
                                               
                                               
  for(a=0; a<coord.n; a++)  if(coord.is_emitter(a)) {
    if(spin_==0)  gg.alloc(lmax_eff,1);        
    else          gg.alloc(lmax_eff,2);        
    ekRa= 4*pi *
          exp(  i_c * (  k0x*(coord.x[a]-coord.x0)
                       + k0y*(coord.y[a]-coord.y0)
                       + k0z*(coord.z[a]-coord.z0))
               - (coord.zs-coord.z[a])*imag(scat.k)/ABS(cos(th0))); 
    for(l=lm=0; l<=lmax_eff; l++)
    for(m=-l; m<=l; m++, lm++) {
      psi=ekRa*conj(Ylm[lm]);
      if(spin_==1 || spin_==0) {                   
          gg.re[lm]    =real(psi);
          gg.im[lm]    =imag(psi);
      } else {
          gg.re[lm+dim]=real(psi);
          gg.im[lm+dim]=imag(psi);
    } }
    scat.R_1TR(gg, a);                
    iter.o.g[a].copy(gg);  gg.free();
  }  else  if(spin_==0) iter.o.g[a].alloc(-1,1);
           else         iter.o.g[a].alloc(-1,2);
  delete [] Ylm;
  for(a=0; a<coord.n; a++)  iter.n.g[a].copy(iter.o.g[a]);
  return 1;
}
void electron_bin::surface_lattice(numero a, numero b, numero al, numero be)
{
  leed_a=a;      leed_b=b;                  
  ax=a*cos(al);  bx=b*cos(al+be);           
  ay=a*sin(al);  by=b*sin(al+be);           
  numero a_b=ax*by-ay*bx;                   
  Gax= 2*pi * by/a_b;
  Gay=-2*pi * bx/a_b;
  Gbx=-2*pi * ay/a_b;
  Gby= 2*pi * ax/a_b;
  area=ABS(a_b);
}
void electron_bin::init_leed_emitters(void)
{
  int i;
  numero ra,rb;                             
  numero ab=ax*bx+ay*by;                    
  numero c1,c2;                             
  numero den=sqr(leed_a*leed_b)-sqr(ab);
  coord.free_emitters();                    
  coord.sort(-1);                           
  for(i=0; i<coord.n; i++) {                
                                            
    ra= (coord.x[i]-leed_x)*ax              
       +(coord.y[i]-leed_y)*ay;             
    rb= (coord.x[i]-leed_x)*bx              
       +(coord.y[i]-leed_y)*by;
    c1=leed_a*(ra*sqr(leed_b)-rb*ab)/den;   
    c2=leed_b*(rb*sqr(leed_a)-ra*ab)/den;   
    if(-coord.eps_d<c1 && c1<leed_a-coord.eps_d &&
       -coord.eps_d<c2 && c2<leed_b-coord.eps_d)  coord.add_emitter(i);
} }
numero electron_bin::intensity(complex *wf)
{
  numero val=0;
  int c;
  if(spin_flag)  val=mod2(wf[0]);
  else  for(c=0; c<n_c; c++) val+=mod2(wf[c]);
  return val;
}
void electron_bin::wf_far(int i, int j, complex *wf)
{
  complex *Ylm;  Ylm=final.Ylm(i,j);
  wf_far(final.th[i], final.fi[j], Ylm, wf);
  delete [] Ylm;
}
void electron_bin::wf_far(numero th, numero fi, complex *Ylm, complex *wf)
{
  int a, c;
  numero kfx,kfy,kfz, ilam, kr=real(scat.k);
  complex *wft, ekRa;
  kfx= kr*sin(th)*cos(fi);
  kfy= kr*sin(th)*sin(fi);
  kfz= kr*cos(th);
  ilam=imag(scat.k)/ABS(cos(th));           
  wft=new complex [n_c];
  for(c=0; c<n_c; c++)  wf[c]=0;
  for (a=0; a<coord.n; a++) {
    ekRa= exp( -i_c * (  kfx*(coord.x[a]-coord.x0)
                       + kfy*(coord.y[a]-coord.y0)
                       + kfz*(coord.z[a]-coord.z0))
               -(coord.zs-coord.z[a])*ilam*(1-molecule) );
    wf_far(iter.n.g[a], Ylm, wft);
    for(c=0; c<n_c; c++)  wf[c]+=wft[c]*ekRa;
  }
  for(c=0; c<n_c; c++)  wf[c]=wf[c]/kr;     
  delete [] wft;
}
void electron_bin::wf_far(vector &a, complex *Ylm, complex *wf)
{
  int l,m,lm, c;  numero *ar,*ai;                
                                                 
  if((th_spin!=0 || fi_spin!=0) && a.n_c==2)     
    a.spin_rotation(0, th_spin, pi-fi_spin);     
                                                 
  for(c=0; c<a.n_c; c++) {
    ar=a.re+(c*lp12[a.lmax]);
    ai=a.im+(c*lp12[a.lmax]);   wf[c]=0;
    for(l=lm=0; l<=a.lmax; l++)
    for(m=-l; m<=l; m++, ar++, ai++, lm++) {
      wf[c]+=complex((*ar),(*ai)) * Ylm[lm];
  } }
  if((th_spin!=0 || fi_spin!=0) && a.n_c==2)     
    a.spin_rotation(fi_spin-pi, -th_spin, 0);    
}                                                
void electron_bin::wf_near(int i, int j, complex *wf)
{
  int a, c;                                      
  complex *wft;
  wft=new complex [n_c];
  for(c=0; c<n_c; c++)  wf[c]=0;
  for (a=0; a<coord.n; a++) {
  }
  delete [] wft;
}
void electron_bin::wf_LEED(complex *wf)
{
  int c;
  numero kf;                                            
  numero fi;                                            
  numero kfx= k0x + leed_i*Gax + leed_j*Gbx;            
  numero kfy= k0y + leed_i*Gay + leed_j*Gby;            
  numero kfz2=sqr(real(scat.k))-sqr(kfx)-sqr(kfy);      
  if(kfz2<=0)  {for(c=0; c<n_c; c++) wf[c]=0;  leed_th=0;}
  else {
    cartesian_to_spherical(kfx,kfy,sqrt(kfz2), kf,leed_th,fi);
    complex *Ylm;   Ylm=new complex [lp12[lmax_eff]];   
    final.init_Ylm(leed_th,fi,Ylm);
    wf_far(leed_th, fi, Ylm, wf);                       
    delete [] Ylm;
} }
void electron_bin::init_wf_PD(int jm, int mjms, complex *eps)
{
  int i;
  for(i=0; i<coord.n; i++)  if(i!=emitter)  iter.o.g[i].alloc(-1,n_c);
  if(manual_flag)  init_wf_manual();
  else {
    if(rmat==NULL && CR_flag) {                           
      scat.get_scatterer(coord.atom[emitter]);
      rrm[0]=rmat_el_CR(l0+1);   
      rrm[1]=rmat_el_CR(l0-1);   
    } else
    if(rmat==NULL) {
      scat.get_scatterer(coord.atom[emitter]);
      rrm[0]=rmat_el_MT(l0+1);   
      rrm[1]=rmat_el_MT(l0-1);   
    } else
    if(excitation==0) {
      rrm[0]=rmat[0].val(real(scat.k));
      rrm[1]=rmat[1].val(real(scat.k));  if(l0==0) rrm[1]=0;
    } else
    if(excitation==1) {
      rrm[0]=rmat[0].val(real(scat.k));
      rrm[1]=rmat[1].val(real(scat.k));
      rrm[2]=rmat[2].val(real(scat.k));
      rrm[3]=rmat[3].val(real(scat.k));  if(l0==0) rrm[1]=rrm[3]=0;
    }
    if(ljmj)  init_wf_ljmj(jm,mjms,eps);
    else      init_wf_lmms(jm,mjms,eps);
  }
  for(i=0; i<coord.n; i++)                       
    if(i!=emitter)  iter.n.g[i].alloc(-1,n_c);
    else            iter.n.g[i].copy(iter.o.g[i]);
}
void electron_bin::init_rmat_MT(void)
{
  numero E=-sqr(scat.a->atZ+0.0)/sqr(n0), cc=c_au;  c_au*=1000;   
  free_rmat_MT();
  Rho_dbound(scat.a->atZ, P0,Q0, E, n0,-l0-1);  c_au=cc;
}
complex electron_bin::rmat_el_MT(int l)
{
  if(l<0)  return 0;
  Salvat Dirac;                                
  spline integral;
  numero r, f,g, rr,ps, cc=c_au;  c_au*=1000;  
  int i;
  Dirac.vint(scat.a->rad, scat.a->re, scat.a->n_rad, scat.a->atZ);
  Dirac.dfree(sqr(real(scat.k))/2, ps, -l-1);
  integral.alloc(P0.n);
  for(i=0; i<P0.n; i++) {
    r=P0.x[i];
    Dirac.wf(r,f,g);
    integral.put(i, r, P0.y[i]*f);
  }
  integral.init();
  rr=integral.integ(P0.a, P0.b, 1);          
if(0) {
  printf("--> %g %g %g %d\n", integral.integ(P0.a, P0.b, 1),
                              integral.integ(P0.a, P0.b, 2), ps, l);
  Dirac.dfree(sqr(real(scat.k))/2, ps, -2-1);
  integral.alloc(P0.n);
  for(i=0; i<P0.n; i++) {
    r=P0.x[i];
    Dirac.wf(r,f,g);
    integral.put(i, r, P0.y[i]*f);
  }
  integral.init();
  printf("--> %g %g %g %d\n", integral.integ(P0.a, P0.b, 1),
                              integral.integ(P0.a, P0.b, 2), ps, l+1);
}
  Dirac.free();  integral.free();  c_au=cc;
  return  euler(rr,ps);
}
numero electron_bin::CR_wf(numero r)
{
  numero val=0, red;
  int n;
  if(l0==0)  red=0.75;  else
  if(l0==1)  red=0.7;   else
             red=1;
  for(n=0; n<CR_nterms; n++)
    val+=CR_C[n]/sqrt(fact[2*CR_n[n]])
         * pow(r,CR_n[n]-1.0)
         * pow(2*CR_xhi[n],CR_n[n]+0.5)
         * exp(-CR_xhi[n]*r/red);
  return val;
}
complex electron_bin::rmat_el_CR(int l)
{
  if(l<0)  return 0;
  Salvat Dirac;                                
  spline integral;
  numero r, f,g, rr,ps, cc=c_au;  c_au*=1000;  
  int i, nr=100;
  numero  rmax=10;
  Dirac.vint(scat.a->rad, scat.a->re, scat.a->n_rad, scat.a->atZ);
  Dirac.dfree(sqr(real(scat.k))/2, ps, -l-1);
  integral.alloc(nr);
  for(i=0; i<nr; i++) {
    r=(i+0.5)*rmax/nr;
    Dirac.wf(r,f,g);
    integral.put(i, r, CR_wf(r)*f);    
  }                                    
  integral.init();                     
                                       
  rr=integral.integ(integral.a, integral.b, 2);   
                                                  
  Dirac.free();  integral.free();  c_au=cc;
  return  euler(rr,ps);
}
complex electron_bin::ang_mat_el(int m0_, int l, int m, complex *eps)
{
  numero  gg[3];               
  complex  val;                
                               
  if(m0_-1==m)  gg[0]=Gaunt(l0,1,l,m0_,-1,m);  else  gg[0]=0;
  if(m0_  ==m)  gg[1]=Gaunt(l0,1,l,m0_, 0,m);  else  gg[1]=0;
  if(m0_+1==m)  gg[2]=Gaunt(l0,1,l,m0_, 1,m);  else  gg[2]=0;
  return  sqrt(4*pi/3) * (gg[0]*eps[0] + gg[1]*eps[1] + gg[2]*eps[2]);
}
void electron_bin::init_wf_lmms(int m0_, int ms0_, complex *eps)
{
  int i,l,m,c;                               
  complex val, rm;                           
  iter.o.g[emitter].alloc(l0+1,n_c);
  for(c=i=0; c<n_c; c++)
  for(l=0; l<=iter.o.g[emitter].lmax; l++)
  for(m=-l; m<=l; m++, i++) {
    iter.o.g[emitter].re[i]=iter.o.g[emitter].im[i]=0;
    if((l==l0+1) || (l==l0-1))
    if(excitation==1 || n_c==1 || (ms0_==1 && c==0) || (ms0_==-1 && c==1)) {
      if(l>l0)  if(excitation==1)  rm=rrm[2*c];    else  rm=rrm[0];
      else      if(excitation==1)  rm=rrm[2*c+1];  else  rm=rrm[1];
      rm=real(scat.k)*sqrt(4*pi) * i_l(-l) * rm;    
      if(n_c==1)  rm=rm*sqrt(2);                    
      val=rm*ang_mat_el(m0_,l,m, eps);
      iter.o.g[emitter].re[i]= real(val);
      iter.o.g[emitter].im[i]= imag(val);
} } }
void electron_bin::init_wf_ljmj(int j0_, int mj0_, complex *eps)
{
  int i,l,m,m0_, dim;                        
  numero cgu = CGspin1_2(2*l0,j0_, 1,mj0_),  
         cgd = CGspin1_2(2*l0,j0_,-1,mj0_);  
  complex val, rm;                           
                                             
  if(2*l0-1!=j0_ && 2*l0+1 != j0_)
    on_error(foutput,"scan pd","unphysical (l0,j0) combination");
  iter.o.g[emitter].alloc(l0+1,n_c);
  dim=lp12[iter.o.g[emitter].lmax];
  for(l=i=0; l<=iter.o.g[emitter].lmax; l++)
  for(m=-l; m<=l; m++, i++) {
    iter.o.g[emitter].re[i]=iter.o.g[emitter].im[i]=0;
    iter.o.g[emitter].re[i+dim]=iter.o.g[emitter].im[i+dim]=0;
    if((l==l0+1)||(l==l0-1)) {
      if(l>l0)  rm=rrm[0];  else  rm=rrm[1];
      rm=real(scat.k)*sqrt(4*pi) * i_l(-l) * rm;  
      m0_=(mj0_-1)/2;                                
      val=rm*ang_mat_el(m0_,l,m, eps);
      iter.o.g[emitter].re[i]    = cgu * real(val);
      iter.o.g[emitter].im[i]    = cgu * imag(val);
      m0_=(mj0_+1)/2;                                
      val=rm*ang_mat_el(m0_,l,m, eps);
      iter.o.g[emitter].re[i+dim]= cgd * real(val);
      iter.o.g[emitter].im[i+dim]= cgd * imag(val);
} } }
void electron_bin::init_wf_manual(void)
{                                          
  int i,c,l,m;
  vector gg;
  iter.o.g[emitter].alloc(lmax_manual,n_c);
  for(c=i=0; c<n_c; c++)
  for(l=0; l<=iter.o.g[emitter].lmax; l++)
  for(m=-l; m<=l; m++, i++) {
    iter.o.g[emitter].re[i] = real(manual[i]);
    iter.o.g[emitter].im[i] = imag(manual[i]);
} }
void electron_bin::read_rmat(FILE *fin)
{
  FILE *frmat;  int i,j,n_k, f;  numero k, x,y,R;
  char units1[str_length_units];        
  char index2[str_length_comments];
  char my_file[str_length_files];
  free_rmat();
  read_name(fin,my_file);
  if(!strcmpC(my_file,"inline"))  frmat=fin;
  else                            frmat=open_file(foutput,my_file,"r");
  n_k=read_int(frmat);
  if(n_k>1)  read_name(frmat,units1);
  read_name(frmat,index2);
  
  
  
  
  
  if(!strcmpC(index2,"regular1")) {excitation=0;  n_rmat=2;  f=1;}  else
  if(!strcmpC(index2,"regular2")) {excitation=0;  n_rmat=2;  f=0;}  else
  if(!strcmpC(index2,"ex1"))      {excitation=1;  n_rmat=4;  f=1;}  else
  if(!strcmpC(index2,"ex2"))      {excitation=1;  n_rmat=4;  f=0;}  else
  on_error(foutput,"read_rmat","wrong argument:",index2);
  if(excitation==1)  ljmj=0;
  l0=read_int(frmat);
  rmat=new splinec [n_rmat];
  for(i=0; i<n_rmat; i++)
    {rmat[i].n=0; rmat[i].alloc(n_k,interpolation_type);}
  for(j=0; j<n_k; j++) {
    if(n_k>1) {k=read_numero(frmat);  k=units_k(units1,k);}
    for(i=0; i<n_rmat; i++) {
      x=read_numero(frmat);  y=read_numero(frmat);
      if(f) {R=x;  x=R*cos(y);  y=R*sin(y);}
      if(n_k>1) rmat[i].put(j,k,complex(x,y));
      else      rmat[i].init(complex(x,y),0);
  } }
  if(n_k>1)  for(i=0; i<n_rmat; i++)  rmat[i].init();
  if(strcmpC(my_file,"inline"))  fclose(frmat);
}
void electron_bin::free_rmat(void)
{
  int i;
  if(rmat!=NULL)  {
    for(i=0; i<n_rmat; i++) rmat[i].free();
    delete [] rmat;  rmat=NULL;
} }
int electron_bin::scan_PD(char *name)
{
  int k,n, emit_s, initial,pol, sct,emit_MT;
  int ljmj_old=ljmj, ms0_old=ms0;
  
  if(manual_flag==1 && (calc.mode_sample>0 || calc.mode_beta>0))
    on_error(foutput,"scan pd",
             "cannot use wf manual with this type of mobility");
  if(manual_flag==1 && manual==NULL)
    on_error(foutput,"scan pd", "empty list wf manual coefficients");
  
  if(manual_flag)                    n_c=n_c_manual;  else
  if(scat.get_n_c()==2)              n_c=2;           else
  if(rmat!=NULL && excitation==1)    n_c=2;           else
  if(ljmj && (j0<1000 || mj0<1000))  n_c=2;           else
  n_c=1;
  if(n_c==1)  {ms0=1;  ljmj=0;}
  scat.set_lmax();
  
  int  jmmin,jmmax,jmstep,jmcomp, mjmsmin,mjmsmax,mjmsstep,mjmscomp,
       jm,                       
       mjms;                     
                                 
  if(manual_flag) {              
    jmmin=1;       jmmax=1;       jmstep=1;    jmcomp=1;
    mjmsmin=1;     mjmsmax=1;     mjmsstep=1;  mjmscomp=1;
  }  else
  if(ljmj==0) {             
    jmmin=-l0;     jmmax=l0;      jmstep=1;    jmcomp=m0;
    mjmsmin=-1;    mjmsmax=1;     mjmsstep=2;  mjmscomp=ms0;
  } else {                  
    jmmin=2*l0-1;  jmmax=2*l0+1;  jmstep=2;    jmcomp=j0;
    mjmsmin=-j0;   mjmsmax=j0;    mjmsstep=2;  mjmscomp=mj0;
  }
  
  calc.n_c=n_c;  calc.select_wf=0;
  calc.init_final_dist();
  calc.init(name, "scan PD");
  
  for(k=0; k<calc.n_k; k++) {    
    scat.set_k(k);               
    final.init_Ylm();
    calc.alloc(coord.n_emit*((jmmax-jmmin)/jmstep+1)
                           *((mjmsmax-mjmsmin)/mjmsstep+1));
    
    
    for(n=initial=0; n<coord.n_emit; n++) {               
      coord.sort(coord.emit[n]);                          
      emitter=coord.emit[n];                              
      emit_s=coord.atom[emitter];                         
      scat.get_scatterer(emit_s);
      if(manual_flag==0 && rmat==NULL) {                  
        if(scat.a->select!=3)
          on_error(foutput,"scan pd", "undefined radial mat. elem.");
        if(CR_flag==0)  init_rmat_MT();                   
      }
      if(scat.a->select==3) {                             
        emit_MT=1;                                        
        sct=1000;  while(scat.get_scatterer(sct,0)==0) sct++;
        coord.atom[emitter]=sct;
        scat.add_screened_scatterer(sct, emit_s);
      }  else  {emit_MT=0;  coord.atom[emitter]=coord.emit_s[n];}
      iter.alloc();                                       
      if(iter.method==3)  scat.get_exact(n_c,lmax_eff);   
      
      for(jm=jmmin; jm<=jmmax; jm+=jmstep)                
      if(jmcomp>1000 || jm==jmcomp) {
        for(mjms=mjmsmin; mjms<=mjmsmax; mjms+=mjmsstep)  
        if(mjmscomp>1000 || mjms==mjmscomp)
        
        for(pol=0; pol<calc.mode_eps; pol++, initial++) { 
            if(manual_flag==0)  calc.set_polarization(pol);
            if(calc.wf_fin==NULL) init_wf_PD(jm, mjms, calc.eps);
            calc.evaluation(initial);
            iter.clear();
      } }                        
      if(iter.method==3)  scat.free_exact();              
      iter.free();
      coord.atom[emitter]=emit_s;
      if(emit_MT) scat.free_scatterer(sct);
    }                            
    scat.free_Green_functions();
    scat.free_t_matrix();
    calc.write(k);               
    final.free_Ylm();
  }                              
  calc.free();                   
  scat.free_rotations();
  free_lmax();
  free_rmat_MT();  ljmj=ljmj_old;  ms0=ms0_old;
  return 0;
}
int electron_bin::scan_MS(char *name, int nwf)
{
  numero *kk;  kk=copy(calc.kk,calc.n_k);  delete [] calc.kk;
  int iwf, n_k=calc.n_k;  calc.n_k=1;  calc.kk=new numero [1];
  char k_units[str_length_units];
  strcpy(k_units, calc.k_units);  strcpy(calc.k_units, "E(eV)");
  
  calc.n_c=wf_ns;  calc.select_wf=0;  calc.ms=1;
  calc.init_final_dist();
  calc.init(name, "scan MS");
  
  for(iwf=0; iwf<nwf; iwf++) {   
    calc.alloc(1);               
    
    iter.alloc();                                         
    if(iter.method==3)  scat.get_exact(wf_ns,lmax_eff);   
    calc.evaluation(0);
    iter.clear();
    if(iter.method==3)  scat.free_exact();                
    iter.free();
    scat.free_Green_functions();
    scat.free_t_matrix();
    calc.write(0);
    final.free_Ylm();
    scat.free_rotations();
    free_lmax();
  }
  calc.free();
  coord.free();
  calc.kk=copy(kk,n_k);  delete [] kk;  calc.n_k=n_k;  calc.ms=0;
  strcpy(calc.k_units, k_units);  strcpy(calc.wf_in,"none");
  return 0;
}
int electron_bin::scan_DLEED(char *name)
{
  int k;
  n_c=1;                    
  scat.set_lmax();               
  coord.sort(-1);
  
  calc.n_c=n_c;  calc.select_wf=0;
  calc.init_final_dist();
  calc.init(name, "scan DLEED");
  iter.alloc();
  
  for(k=0; k<calc.n_k; k++) {    
    scat.set_k(k);               
    final.init_transmission(0);  
                                 
    final.init_Ylm();
    calc.alloc(1);               
    if(calc.wf_fin==NULL) init_wf_LEED(0);
    if(iter.method==3)  scat.get_exact(n_c, lmax_eff);  
    calc.evaluation(0);
    if(iter.method==3)  scat.free_exact();              
    iter.clear();
    scat.free_Green_functions();
    scat.free_t_matrix();
    calc.write(k);               
    final.free_Ylm();
  }
  calc.free();                   
  iter.free();
  scat.free_rotations();
  free_lmax();
  return 0;
}
int electron_bin::scan_LEED(char *name)
{
  int k;
  n_c=1;                    
  scat.set_lmax();               
  coord.sort(-1);
  
  calc.n_c=n_c;  calc.select_wf=-2;   
  calc.init_final_f(1);               
  calc.init(name, "scan LEED");
  iter.alloc();
  
  for(k=0; k<calc.n_k; k++) {    
    scat.set_k(k);               
    final.init_transmission(0);  
                                 
    calc.alloc(1);               
    if(calc.wf_fin==NULL) init_wf_LEED(0);
    if(iter.method==3)  scat.get_exact(n_c,lmax_eff);  
    calc.evaluation(0);
    if(iter.method==3)  scat.free_exact();             
    iter.clear();
    scat.free_Green_functions();
    scat.free_t_matrix();
    calc.write(k);               
  }
  calc.free();                   
  iter.free();
  scat.free_rotations();
  free_lmax();
  return 0;
}
void electron_bin::program(char *my_file)
{
  FILE *fprog=NULL, *ftmp=NULL;           
  char command[str_length_comments];      
  char light[str_length_comments];        
  char units[str_length_units];           
  char name[str_length_files];            
  numero  al,be,del, a1,a2, xx,yy,zz, th,fi,  tt,ff, val1,val2;
  int     *o, e1,e2, i,j,k,l,m,n;
  complex *kc;
  char    c;
  numero  x0,y0,z0, x1,y1,z1, x2,y2,z2, r, th1,th2, fi1,fi2;
  int     n_1, n_2;
  if(!strcmpC(my_file,"inline"))  fprog=stdin;  else
  if((fprog=fopen(my_file,"r"))==NULL)
      fprog=open_file(foutput,strcat(my_file,".ms"),"r");
  for(read_name(fprog, command);
      strcmpC(command,"end") &&
      strcmpC(command,"quit") &&
      strcmpC(command,"exit");
      read_name(fprog, command)) {
    if(command[0]=='"') {
      i=strlen(command);
      if(i>1) {
        for(j=1; j<i-1; j++)  {
          c=command[j];
          if(c=='\\') {
            c=command[++j];
            if(c=='\\' || c=='"')  fprintf(foutput,"%c",c);  else
            if(c=='n')  fprintf(foutput,"\n");
          }  else  fprintf(foutput,"%c",c);
        }
        if(j==i)  fscanf(fprog,"%c",&c);  else  c=command[i-1];
      }  else  fscanf(fprog,"%c",&c);
      while (c!='"') {
        if(c=='\\') {
          fscanf(fprog,"%c",&c);
          if(c=='\\' || c=='"')  fprintf(foutput,"%c",c);  else
          if(c=='n')  fprintf(foutput,"\n");
        }  else  fprintf(foutput,"%c",c);
        fscanf(fprog,"%c",&c);
    } }  else
    if((i=algebra_locate(command,'='))>0) {   
      command[i]=0;
      algebra_define(command,command+i+1);
    } else
    if(!strcmpC(command,"add")) {
      read_name(fprog,name);
      if(!strcmpC(name,"cluster")) {
        read_name(fprog,name);                
        i=coord.read_atomic_type(fprog);
        read_name(fprog,units);  val1=units_conversion_factor(units);
        xx=val1*alread_numero(fprog);
        yy=val1*alread_numero(fprog);
        zz=val1*alread_numero(fprog);
        if(!strcmpC(name,"atom")) {           
          if(coord.oriented) {
            th=pi/180*alread_numero(fprog);  fi=pi/180*alread_numero(fprog);}
          else  th=fi=0;
          coord.add_atom(xx,yy,zz, th,fi, i);
        }  else
        if(!strcmpC(name,"half-layer")) {     
          a1=val1*alread_numero(fprog);    a2=val1*alread_numero(fprog);
          al=pi/180*alread_numero(fprog);  del=pi/180*alread_numero(fprog);
          ff=pi/180*alread_numero(fprog);
          if(coord.oriented) {
            th=pi/180*alread_numero(fprog);  fi=pi/180*alread_numero(fprog);}
          else  th=fi=0;
          coord.add_half_layer(xx,yy,zz, a1,a2, al,del, ff, th,fi, i);
        }  else
        if(!strcmpC(name,"layer")) {          
          a1=val1*alread_numero(fprog);      a2=val1*alread_numero(fprog);
          al=pi/180*alread_numero(fprog);    del=pi/180*alread_numero(fprog);
          if(coord.oriented) {
            th=pi/180*alread_numero(fprog);  fi=pi/180*alread_numero(fprog);}
          else  th=fi=0;
          coord.add_layer(xx,yy,zz, a1,a2, al,del, th,fi, i);
        }  else
        if(!strcmpC(name,"row")) {            
          a1=val1*alread_numero(fprog);
          tt=pi/180*alread_numero(fprog);    ff=pi/180*alread_numero(fprog);
          if(coord.oriented) {
            th=pi/180*alread_numero(fprog);  fi=pi/180*alread_numero(fprog);}
          else  th=fi=0;
          coord.add_row(xx,yy,zz, a1, tt,ff, th,fi, i);
        }  else
        if(!strcmpC(name,"surface")) {        
          a1=val1*alread_numero(fprog);
          read_name(fprog,name);
          if(!strcmpC(name,"hcp0001")) a2=val1*alread_numero(fprog); else a2=0;
          if(coord.oriented) {
            th=pi/180*alread_numero(fprog);  fi=pi/180*alread_numero(fprog);}
          else  th=fi=0;
          coord.add_surface(xx,yy,zz, a1,a2, name, th,fi, i);
        }  else
        on_warning(foutput,"(input) add cluster","command ignored:", name);
      } else  on_error(foutput,"(input) add","wrong keyword:", name);
    }  else
    if(!strcmpC(command,"attenuation-type"))
      scat.attenuation_type=alread_int(fprog);  else
    if(!strcmpC(command,"beamline"))  {
      read_name(fprog,command);
      if(!strcmpC(command,"ALS-7.0.1"))  i=1;      
      if(!strcmpC(command,"Aebi"))       i=2;      
      if(i==2) {
          calc.beamline_th=pi/180*alread_numero(fprog);   
          calc.beamline_fi=pi/180*alread_numero(fprog);   
      }
      calc.init_beamline(i);
    }  else
    if(!strcmpC(command,"beta")) calc.beta=pi/180*alread_numero(fprog); else
    if(!strcmpC(command,"clear")) {
      read_name(fprog,command);
      if(!strcmpC(command,"cluster"))                  
        {coord.free();  scat.free_rotations();}   else
      if(!strcmpC(command,"emitters"))                 
        coord.free_emitters();                    else
      if(!strcmpC(command,"imfp"))                     
        {scat.iimfp.free();  scat.iimfp.n=0;}     else
      if(!strcmpC(command,"rmat"))                     
        free_rmat();                              else
      if(!strcmpC(command,"scatterer"))                
        scat.free_scatterer(alread_int(fprog));   else
      if(!strcmpC(command,"scatterers"))               
        scat.free_scatterers();                   else
        on_warning(foutput,"(input) clear","command ignored:", command);
    }  else
    if(!strcmpC(command,"cluster")) {
      read_name(fprog,name);
      if(!strcmpC(name,"input")) {            
        read_name(fprog,name);
        coord.read(fprog,name);               
      } else                                  
      if(!strcmpC(name,"output")) {
        read_name(fprog,units);  read_name(fprog,name);
        coord.sort(-1);
        coord.write(name,units);
      } else
      if(!strcmpC(name,"natoms")) {           
        coord.nnn=alread_int(fprog);  coord.set_n();
      } else                                  
      if(!strcmpC(name,"reference-point")) {
        read_name(fprog,units);
        coord.x0=units_conversion_factor(units)*alread_numero(fprog);
        coord.y0=units_conversion_factor(units)*alread_numero(fprog);
        coord.z0=units_conversion_factor(units)*alread_numero(fprog);
      }  else                                 
      if(!strcmpC(name,"Rmax")) {
        read_name(fprog,units);
        coord.Rmax=units_conversion_factor(units)*alread_numero(fprog);
      } else                                  
      if(!strcmpC(name,"surface")) {          
        read_name(fprog,units);
        if(!strcmpC(units,"on"))   coord.surface=1;  else
        if(!strcmpC(units,"off"))  coord.surface=0;  else
        coord.zs=units_conversion_factor(units)*alread_numero(fprog);
      } else  on_error(foutput,"(input) cluster","wrong keyword:", name);
    }  else
    if(!strcmpC(command,"Debye-temperature")) {
      read_name(fprog,name);                         
      l=element_properties(name);
      if(l==0)  l=alread_int(name);
      val1=alread_numero(fprog);  if(val1<1)  val1=1;  
      scat.T_Debye[l]=val1;
    }  else
    if(!strcmpC(command,"displace")) {
      read_name(fprog,name);
      if(!strcmpC(name,"cluster")) {         
        read_name(fprog,units);
        xx=units_conversion_factor(units)*alread_numero(fprog);
        yy=units_conversion_factor(units)*alread_numero(fprog);
        zz=units_conversion_factor(units)*alread_numero(fprog);
        coord.displace(xx,yy,zz);
      } else
      on_warning(foutput,"(input) displace","command ignored:", name);
    }  else
    if(!strcmpC(command,"dmax")) {
      read_name(fprog,units);
      scat.dmax=units_conversion_factor(units)*alread_numero(fprog);
    }  else
    if(!strcmpC(command,"emission")) {
      read_name(fprog, command);
      if(!strcmpC(command,"energy")) {
        read_name(fprog, units);       
        val1=alread_numero(fprog);       
        val2=alread_numero(fprog);       
        n=alread_int(fprog);             
        calc.init_k(val1,val2,n,units);
      } else
      if(!strcmpC(command,"angle")) {
        read_name(fprog, command);
        if(!strcmpC(command,"theta")) {
          val1=pi/180*alread_numero(fprog);   
          val2=pi/180*alread_numero(fprog);   
          n=alread_int(fprog);                
          final.init_th(val1,val2,n);       
        } else
        if(!strcmpC(command,"phi")) {
          val1=pi/180*alread_numero(fprog);   
          val2=pi/180*alread_numero(fprog);   
          n=alread_int(fprog);                
          final.init_phi(val1,val2,n);      
        } else
        if(!strcmpC(command,"normal")) {
          final.init_th(0,0,1);  final.init_phi(0,0,1);
        } else
        if(!strcmpC(command,"window")) {
          calc.thave=pi/180*alread_numero(fprog);   
        }
        else  on_error(foutput,"(input) emission angle",
                       "theta/phi/normal expected instead of", command);
      } else  on_error(foutput,"(input) emission",
                       "energy/angle expected instead of", command);
    }  else
    if(!strcmpC(command,"emitters")) {
      read_name(fprog,name);
      if(!strcmpC(name,"all")) {                                
        e1=alread_int(fprog);
        e2=alread_int(fprog);
        coord.all_emitters(e1,e2);
      } else                                                    
      if(!strcmpC(name,"plane")) {
        read_name(fprog,units);
        xx=units_k(units, alread_numero(fprog));
        yy=units_k(units, alread_numero(fprog));
        zz=units_k(units, alread_numero(fprog));
        th=pi/180*alread_numero(fprog);
        fi=pi/180*alread_numero(fprog);
        coord.add_emitter_plane(xx,yy,zz,th,fi);
      } else                                                    
      if(!strcmpC(name,"column")) {
        read_name(fprog,units);
        xx=units_k(units, alread_numero(fprog));
        yy=units_k(units, alread_numero(fprog));
        coord.add_emitter_column(xx,yy);
      }  else {                                                 
        n=alread_int(name);   read_name(fprog,units);
        for(i=0; i<n; i++) {
          if(units[0]=='l') {                                   
            xx=units_k(units, alread_numero(fprog));
            yy=units_k(units, alread_numero(fprog));
            zz=units_k(units, alread_numero(fprog));
            val1=0.05/a0_au;
            if(val1<coord.eps_d)  val1=coord.eps_d;
            e1=coord.get_atom(xx,yy,zz, val1);
            if(e1==-1)  on_error(foutput,"(input) emitters",
                                 "cluster position not found");
          }  else  if(i==0)  e1=alread_int(units)-1;
                   else      e1=alread_int(fprog)-1;              
          e2=alread_int(fprog);
          coord.add_emitter(e1,e2);
    } } }  else
    if(!strcmpC(command,"fixed")) {
      read_name(fprog,name);             
      calc.init_mode_sample(0);          
      calc.init_mode_beta(0);            
                                         
    } else                               
    if(!strcmpC(command,"giant-matrix-output")) {
      read_name(fprog,scat.aa_output);
      scat.aa_output_flag=1;
      iter.method= 3;
    }  else
    if(!strcmpC(command,"imfp")) {
      scat.iimfp.free();  scat.iimfp.n=0;
      read_name(fprog,name);
      if(!strcmpC(name,"TPP-2M")) {          
        scat.TPP_rho=alread_numero(fprog);     
        scat.TPP_Nv=alread_numero(fprog);      
        scat.TPP_Ep=alread_numero(fprog);      
        scat.TPP_Eg=alread_numero(fprog);      
        if(scat.TPP_rho<=0 || scat.TPP_Nv<=0
                           || scat.TPP_Ep<=0 || scat.TPP_Eg<0)
          on_error(foutput,"(input) imfp TPP-2M", "wrong parameters");
        scat.iimfp_flag=1;
      } else {
        scat.read_imfp(fprog,name);
        scat.iimfp_flag=0;
    } }  else
    if(!strcmpC(command,"incidence")) {
      read_name(fprog, light);
      if(!strcmpC(light,"normal"))  calc.init_incidence(0,0);
      else {
        val1=pi/180*alread_numero(light);
        val2=pi/180*alread_numero(fprog);
        calc.init_incidence(val1,val2);
    } }  else
    if((!strcmpC(command,"include")) || (!strcmpC(command,"insert")))  {
      read_name(fprog,name);
      program(name);
    }  else
    if(!strcmpC(command,"initial")) {
      read_name(fprog,command);
      CR_flag=0;                                      
      if(!strcmpC(command,"lmms")) ljmj=0;               else
      if(!strcmpC(command,"ljmj")) ljmj=1;               else
      if(!strcmpC(command,"n0"))   n0 =alread_int(fprog);  else
      if(!strcmpC(command,"l0"))   l0 =alread_int(fprog);  else
      if(!strcmpC(command,"m0"))   m0 =alread_int(fprog);  else
      if(!strcmpC(command,"j0"))   j0 =alread_int(fprog);  else
      if(!strcmpC(command,"ms0"))  ms0=alread_int(fprog);  else
      if(!strcmpC(command,"mj0"))  mj0=alread_int(fprog);  else
      if(!strcmpC(command,"state")) {
        read_name(fprog,command);
        if(!strcmpC(command,"C-R")) {                 
          
          CR_nterms=alread_int(fprog);  CR_flag=1;
          for(j=0; j<CR_nterms; j++) {
            CR_n[j]=alread_int(fprog);
            CR_xhi[j]=alread_numero(fprog);
            CR_C[j]=alread_numero(fprog);
	  }
	} else {
          j=element_orbital(command, n,k);
          if(j==-1)
            on_error(foutput,"(input) initial state", "wrong state:",command);
          n0=n;
          if(k<0) l0=-k-1; else l0=k;
          if(j<1000)  if(k<0) j0=2*l0+1; else j0=2*l0-1;
	}
      }  else  on_error(foutput,"(input) initial", "wrong keyword:",command);
    }  else
    if(!strcmpC(command,"interpolation")) {
      read_name(fprog,command);
      if(!strcmpC(command,"spline"))   interpolation_type=0;  else
      if(!strcmpC(command,"linear"))   interpolation_type=1;  else
        on_warning(foutput,"(input) interpolation",
                   "command ignored:",command);
    }  else
    if(!strcmpC(command,"iteration")) {
      read_name(fprog,command);
      if(!strcmpC(command,"recursive"))  iter.method=-1;    else
      if(!strcmpC(command,"recursion"))  iter.method=-1;    else
      if(!strcmpC(command,"Haydock"))    iter.method=-1;    else
      if(!strcmpC(command,"Jacobi"))     iter.method= 0;    else
      if(!strcmpC(command,"Gauss"))      iter.method= 1;    else
      if(!strcmpC(command,"SR"))         iter.method= 2;    else
      if(!strcmpC(command,"exact"))      iter.method= 3;    else
      on_warning(foutput,"(input) iteration","command ignored:",command);
      if(!strcmpC(command,"SR"))         iter.eta=alread_numero(fprog);
    }  else
    if(!strcmpC(command,"leed")) {
      read_name(fprog,command);
      if(!strcmpC(command,"beam")) {     
        leed_i=alread_int(fprog);
        leed_j=alread_int(fprog);
      }  else
      if(!strcmpC(command,"emitters"))   
        init_leed_emitters();
      else
      if(!strcmpC(command,"lattice")) {  
        read_name(fprog,units);
        a1=units_k(units, alread_numero(fprog));
        a2=units_k(units, alread_numero(fprog));
        al=pi/180*alread_numero(fprog);
        be=pi/180*alread_numero(fprog);
        surface_lattice(a1,a2,al,be);
      }  else
      if(!strcmpC(command,"origin")) {   
        read_name(fprog,units);
        leed_x=units_k(units, alread_numero(fprog));
        leed_y=units_k(units, alread_numero(fprog));
      } }  else
    if(!strcmpC(command,"lmax"))  lmax_max=alread_int(fprog);  else
    if(!strcmpC(command,"mesh1")) {
      read_name(fprog,units);
      x0=units_k(units,alread_numero(fprog));
      y0=units_k(units,alread_numero(fprog));
      z0=units_k(units,alread_numero(fprog));
      x1=units_k(units,alread_numero(fprog));
      y1=units_k(units,alread_numero(fprog));
      z1=units_k(units,alread_numero(fprog));
      x2=units_k(units,alread_numero(fprog));
      y2=units_k(units,alread_numero(fprog));
      z2=units_k(units,alread_numero(fprog));
      n_1=alread_int(fprog);
      n_2=alread_int(fprog);
      final.init_mesh(x0,y0,z0, x1,y1,z1, x2,y2,z2, n_1,n_2);
    }  else
    if(!strcmpC(command,"mesh2")) {
      read_name(fprog,units);
      r=units_k(units,alread_numero(fprog));
      th1=pi/180*alread_numero(fprog);
      th2=pi/180*alread_numero(fprog);
      fi1=pi/180*alread_numero(fprog);
      fi2=pi/180*alread_numero(fprog);
      n_1=alread_int(fprog);
      n_2=alread_int(fprog);
      final.init_mesh(r, th1,th2, fi1,fi2, n_1,n_2);
    }  else
    if(!strcmpC(command,"momenta")) {
      read_name(fprog, units);         
      n=alread_int(fprog);               
      kc=new complex [n];
      for(i=0; i<n; i++) {
        val1=units_k(units,alread_numero(fprog));
        val2=units_k(units,alread_numero(fprog));
        kc[i]=complex(val1,val2);
      }
      calc.init_k(kc,n,units);
    }  else
    if(!strcmpC(command,"molecule")) {
      molecule=1;
    }  else
    if(!strcmpC(command,"movable")) {
      read_name(fprog,name);             
      calc.init_mode_sample(1);          
      calc.init_mode_beta(1);            
                                         
    } else                               
    if(!strcmpC(command,"muffin-tin"))  scat.set_muffin_tin_potentials();
    else
    if(!strcmpC(command,"orders")) {
      n=alread_int(fprog);  o=new int [n];
      for(i=0; i<n; i++)  o[i]=alread_int(fprog);
      calc.init_orders(n,o);
    }  else
    if(!strcmpC(command,"orientation")) {
      read_name(fprog,name);
      if(!strcmpC(name,"on"))   coord.oriented=1;   else
      if(!strcmpC(name,"off"))  coord.oriented=0;  else {
        a1=pi/180*alread_numero(name);
        a2=pi/180*alread_numero(fprog);
        for(i=0; i<coord.nn; i++) {coord.th[i]=a1;  coord.fi[i]=a2;}
    } }  else
    if(!strcmpC(command,"output"))  {
      close_file(foutput);
      read_name(fprog,name);
      foutput=open_file(foutput,name,"a");
    }  else
    if(!strcmpC(command,"periodic")) {
      read_name(fprog,command);
      if(!strcmpC(command,"off"))   scat.separation= 1;  else
      if(!strcmpC(command,"on"))    scat.separation=-1;  else
      if(!strcmpC(command,"test"))  scat.separation= 0;  else
        on_error(foutput,"(input) periodic","wrong keyword:",command);
    }  else
    if(!strcmpC(command,"polarization")) {
      read_name(fprog, light);
      if(!strcmpC(light,"LPx")) calc.init_polarization(0,    0);     else
      if(!strcmpC(light,"LPy")) calc.init_polarization(pi/2, 0);     else
      if(!strcmpC(light,"LCP")) calc.init_polarization(pi/4, pi/2);  else
      if(!strcmpC(light,"RCP")) calc.init_polarization(pi/4,-pi/2);  else
      if(!strcmpC(light,"direction")) {
        val1=pi/2-pi/180*alread_numero(fprog);      
        val2=pi  +pi/180*alread_numero(fprog);      
        normalize_angles(val1,val2);
        calc.init_incidence(val1,val2);
        calc.init_polarization(0,0);
      }  else {
       val1=pi/180*alread_numero(light);
       val2=pi/180*alread_numero(fprog);
       calc.init_polarization(val1,val2);
    } }  else
    if(!strcmpC(command,"remove")) {
      read_name(fprog,name);
      if(!strcmpC(name,"cluster")) {          
        read_name(fprog,name);
        if(!strcmpC(name,"non-emitters"))     
          coord.delete_non_emitters();
        else {
          read_name(fprog,units);  val1=units_conversion_factor(units);
          xx=val1*alread_numero(fprog);  yy=val1*alread_numero(fprog);
          zz=val1*alread_numero(fprog);
          if(!strcmpC(name,"atom"))             
            coord.delete_atom(xx,yy,zz);
          else
          if(!strcmpC(name,"layer")) {          
            a1=val1*alread_numero(fprog);  a2=val1*alread_numero(fprog);
            al=pi/180*alread_numero(fprog);    del=pi/180*alread_numero(fprog);
            coord.delete_layer(xx,yy,zz, a1,a2, al,del);
          }  else
          if(!strcmpC(name,"row")) {            
            a1=val1*alread_numero(fprog);
            tt=pi/180*alread_numero(fprog);    ff=pi/180*alread_numero(fprog);
            coord.delete_row(xx,yy,zz, a1, tt,ff);
          }  else  on_warning(foutput,"(input) remove cluster",
                                      "command ignored:", name);
        }
      } else  on_error(foutput,"(input) remove","wrong keyword:", name);
    }  else
    if(!strcmpC(command,"report")) {
      read_name(fprog,command);
      if(!strcmpC(command,"emitters"))     coord.report_emitters();  else
      if(!strcmpC(command,"links"))        scat.report_links();      else
      if(!strcmpC(command,"lmax"))
        fprintf(foutput,"--- lmax = %d\n", size_of_lmax);            else
      if(!strcmpC(command,"natoms"))
          fprintf(foutput,
          "--- the cluster consists of %d atoms (%d selected, %d used)\n",
          coord.nn, coord.set_n(), calc.n_atoms_last);               else
      if(!strcmpC(command,"scatterer"))
          scat.scan_scatterer(alread_int(fprog),"none",-1);          else
      if(!strcmpC(command,"size"))         report_size();            else
      on_warning(foutput,"(input) report","command ignored:",command);
    }  else
      if(!strcmpC(command,"rmat")) {
        read_rmat(fprog);
        CR_flag=0;                                    
      } else
    if(!strcmpC(command,"rotate")) {
      read_name(fprog,name);
      if(!strcmpC(name,"cluster")) {         
        read_name(fprog,name);
        if(name[0]=='x')  coord.rotate(0, pi/180*alread_numero(fprog));
        if(name[0]=='y')  coord.rotate(1, pi/180*alread_numero(fprog));
        if(name[0]=='z')  coord.rotate(2, pi/180*alread_numero(fprog));
      } else
      on_warning(foutput,"(input) rotate","command ignored:", name);
    }  else
      
    if(!strcmpC(command,"scan")) {
      read_name(fprog,command);
      if(!strcmpC(command,"integrated")) {                
        calc.integ=2;
        read_name(fprog,command);
      }  else
      if(!strcmpC(command,"half-integrated")) {
        calc.integ=1;
        read_name(fprog,command);
      }  else  calc.integ=0;
      i=-1;
      if(!strcmpC(command,"pd"))            i=0;  else    
      if(!strcmpC(command,"ms"))            i=1;  else    
      if(!strcmpC(command,"dleed"))         i=2;  else    
      if(!strcmpC(command,"leed"))          i=3;          
      if(!strcmpC(command,"scatterer")) {   i=5;          
        l=alread_int(fprog);
        read_name(fprog, name);             j=0;
        if(!strcmpC(name,"phase-shifts"))  {j=2;  read_name(fprog, name);}
        if(!strcmpC(name,"potential"))     {j=3;  read_name(fprog, name);}
      }  else
      if(!strcmpC(command,"imfp")) {        i=6;          
        read_name(fprog, name);
      }  else {
        if(i==1)  n=alread_int(fprog);
        read_name(fprog, name);
        if(!strcmpC(name,"amplitude")) {                  
          read_name(fprog,name);
          calc.n_intensities=0;
        }  else  calc.n_intensities=1;
      }
      if(i==0)  scan_PD(name);                    else
      if(i==1)  scan_MS(name, n);                 else
      if(i==2)  scan_DLEED(name);                 else
      if(i==3)  scan_LEED(name);                  else
      if(i==5)  scat.scan_scatterer(l,name,j);    else
      if(i==6)  scat.scan_imfp(name);             else
      on_error(foutput,"(input) scan","wrong keyword",command);
      scat.aa_output_flag=0;
    }  else
    if(!strcmpC(command,"scatterer")) {
      j=alread_int(fprog);
      read_name(fprog,command);
      if(!strcmpC(command,"mass")) {
        scat.get_scatterer(j);
        read_name(fprog,name);
        if(element_properties(name))  scat.a->mass=element_mass*amu_au;
        else                          scat.a->mass=alread_numero(name);
      }  else
      if(!strcmpC(command,"sphere")) {
        scat.add_scatterer(j);
        read_name(fprog,units);
        scat.a->a=units_conversion_factor(units)*alread_numero(fprog);  
        read_name(fprog,units);
        scat.a->V=units_conversion_factor(units)*alread_numero(fprog);  
        scat.a->select=-1;
        scat.a->lmax=infinite_int;
        scat.a->n_c=1;
      }  else  {
        scat.add_scatterer(j);
        scat.a->read(fprog, command);
    } }  else
    if(!strcmpC(command,"scattering-so")) scat.scattering_so=alread_int(fprog);
    else
    if(!strcmpC(command,"screening-length")) {
      read_name(fprog,units);
      scat.screening_length=units_conversion_factor(units)
                            * alread_numero(fprog);
      if(scat.screening_length<=0)  scat.screening_length=1e-10;
    }  else
    if(!strcmpC(command,"semi-movable")) {
        read_name(fprog,name);
        calc.init_mode_sample(1);
        calc.init_mode_beta(0);              
    }  else
      
    if(!strcmpC(command,"spin"))  {
      read_name(fprog,command);
      if(!strcmpC(command,"on"))   spin_flag=1;  else
      if(!strcmpC(command,"off"))  spin_flag=0;  else {
        th_spin=pi/180*alread_numero(command);
        fi_spin=pi/180*alread_numero(fprog);
    } }  else
    if(!strcmpC(command,"temperature")) {
      scat.T=alread_numero(fprog);
      val1=alread_numero(fprog);  if(val1<1)  val1=1;
      for(i=0; i<n_Debye_temperature_elements; i++)
         scat.T_Debye[i]=val1;
    }  else
    if(!strcmpC(command,"temperature-T"))
      scat.T=alread_numero(fprog);
    else
    if(!strcmpC(command,"time"))  time_print(foutput);  else
    if(!strcmpC(command,"tmat")) {
      l=alread_int(fprog);             
      read_name(fprog,name);         
      calc.scan_tmat(name, l);       
                                     
    }  else                          
    if(!strcmpC(command,"tolerance")) {
      read_name(fprog,command);
      if(!strcmpC(command,"angles"))
        coord.eps_th=pi/180*alread_numero(fprog);       else
      if(!strcmpC(command,"distances"))
        coord.eps_d=alread_numero(fprog);               else
      if(!strcmpC(command,"tmat"))
        scat.epsilon=alread_numero(fprog);              else
        on_error(foutput,"(input) precision","wrong keyword",command);
    }  else
    if(!strcmpC(command,"V0"))  {
      read_name(fprog,command);
      scat.V0=alread_numero(fprog) * units_conversion_factor(command);
    }  else
    if(!strcmpC(command,"verbose")) {
      read_name(fprog,name);
        if(!strcmpC(name,"on"))  calc.verbose=1;  else  calc.verbose=0;
    } else
    if(!strcmpC(command,"Vi")) {
      scat.iimfp.free();  scat.iimfp.n=0;  scat.iimfp_flag=2;
      read_name(fprog, units);       
      scat.Vi=alread_numero(fprog) * units_conversion_factor(units);
    } else
    if(!strcmpC(command,"wf")) {
      read_name(fprog,command);
      
      if(!strcmpC(command,"input"))   {
              wf_ns=alread_int(fprog);  read_name(fprog,calc.wf_in);}  else
      if(!strcmpC(command,"output"))    read_name(fprog,calc.wf_out);  else
      if(!strcmpC(command,"manual")) {      
                                            
        read_name(fprog,name);
        if(!strcmpC(name,"on"))  manual_flag=1;  else
        if(!strcmpC(name,"off")) manual_flag=0;  else {
          if(manual!=NULL) {delete [] manual;  manual=NULL;}
          n_c_manual=alread_int(name);
          lmax_manual=alread_int(fprog);
          manual=new complex [n_c_manual*sqr(lmax_manual+1)];
          for(n=i=0; n<n_c_manual; n++)
          for(l=0; l<=lmax_manual; l++)
          for(m=-l; m<=l; m++, i++) {
            a1=alread_numero(fprog);
            a2=alread_numero(fprog);  manual[i]=complex(a1,a2);
      } } }
    } else
    on_warning(foutput,"input file","command ignored:",command);
  }
  close_file(fprog);
  close_file(foutput);  foutput=stdout;
}
void copyright(void)
{
  printf(
  "---------------------------------------------------------------------\n");
  printf("EDAC: Electron Diffraction in Atomic Clusters\n");
  printf("\n");
  printf("Copyright (c) 1997-2000 F. Javier Garcia de Abajo\n");
  printf("Lawrence Berkeley National Laboratory (LBNL)\n");
  printf("\n");
  printf("Please cite the following article ");
  printf("when EDAC is used in published work:\n");
  printf("\n");
  printf("--- Application to the case of PD and LEED:\n");
  printf("    F. J. Garcia de Abajo, M. A. Van Hove, and C. S. Fadley,\n");
  printf("          Phys. Rev. B. 63, 075404 (2001)\n");
  printf(
  "---------------------------------------------------------------------\n");
  printf("\n");
}
#endif  
int main(int argc, char **argv)
{
  int i,j=1;
  particle_type=electrones;
  if(argc<=2)               copyright();  else
  if(strcmp(argv[1],"-t"))  copyright();  else  j++;
  init_fact();                   
  if(argc==1)  {         printf("Inline mode\n\n");
                         electron.program("inline");
  }  else
  for(i=j; i<argc; i++) electron.program(argv[i]); 
  fprintf(foutput,"That's all, folks!\n");
  return 0;
}
