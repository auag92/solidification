
#include "stdio.h"
#include "math.h"
#include <sys/stat.h>

#define MESHX 200
#define MESHX2 (MESHX*MESHX)
#define deltat (0.02)
#define deltax (1.0)
#define deltax2 (deltax*deltax)


#define K (0.5)             /*Partition Coefficient*/
#define G (1.0)             /*Surface Energy*/
#define M (1.0)             /*Mobility*/
#define E (4.0)            /*epsilon - dimensions of length [m]*/
#define tau (1.0)

#define ntimesteps (100000)
#define saveT (1000)

#define deltaMu (0.05)
#define Mu (1.0)

// #define DIRICHLET
#define CONSTANT_X_0 (0.1)
#define CONSTANT_X_1 (0.5)

#define NEUMANN

// #define PERIODIC

// #define SINE_PROFILE
#ifdef SINE_PROFILE
#define wavenumber 4
#define amp 0.1
#endif

#define STEP_PROFILE

double phi_new[MESHX2], phi_old[MESHX2];
double mu_new[MESHX2], mu_old[MESHX2];
double lap_phi[MESHX2], lap_mu[MESHX2];
double conc[MESHX2];

double inv_deltax2 = (1.0/deltax2);

void initialize();
void boundary(double *c);
void write2file ( long t);
void update();
void laplacian(double *f, double *lap);
double h(double phi);
double h1(double phi);
double DW(double phi);
double Cs(double mu);
double Cl(double mu);
double khai(double phi);
void concentration();

void main() {
  long i, j, z, t;
  double p,dp,du;

  initialize();
  boundary(phi_old);
  boundary(mu_old);

  for (t=0; t < ntimesteps; t++) {

    laplacian(phi_old,     lap_phi);
    laplacian(mu_old,     lap_mu);

    for (i=0; i < (MESHX); i++) {
      
      for (j=0; j < (MESHX); j++){
	
	z= i*MESHX + j;
	p = phi_old[z];
	dp = deltat*(2.0*G*E*lap_phi[z] - (G/E)*18.0*(p)*(1.0-p)*(1.0-2.0*p) + (mu_old[z] - Mu)*(K-1)*(mu_old[z])*6*p*(1-p))/(tau*E);
	phi_new[z] = p + dp;
	du = deltat*M*lap_mu[z] - (K-1)*mu_old[z]*6*p*(1-p)*dp;
	mu_new[z] = mu_old[z]  + du/(1+(K-1)*p*p*(3-2*p));
      }
    }

    boundary(phi_new);
    boundary(mu_new);
    concentration();
    update();
    if((t%saveT) ==0) {
      write2file(t);
    }
  }

  printf("Yay! We are done.\n");
}

void update() {
  long i, j, z;
  for (i=0; i < MESHX; i++) {
    for (j=0; j < MESHX; j++){
      z= i*MESHX + j;
      phi_old[z]=phi_new[z];
      mu_old[z]=mu_new[z];
      
    }
  }
}

void laplacian(double *f, double *lap) {
  long i,j,z;
    
  for (i=1; i< MESHX -1; i++)
  {
    for (j=1; j< MESHX -1; j++)
    {
      
      z= i*MESHX + j;
      lap[z] = (f[z-1] + f[z+1] -2.0*f[z])*inv_deltax2 + (f[z+MESHX] + f[z-MESHX] -2.0*f[z])*inv_deltax2 ;

    }

  }
 }

void concentration(){

  double p,u,h;
  int i, j, z;
  for ( i = 0; i < MESHX; i++)
  {
    for ( j = 0; j < MESHX; j++){
      z= i*MESHX + j;
      p = phi_new[z];
      u = mu_new[z];
      h = p*p*(3-2*p);
      conc[z] = u*(1-h) + (h)*K*u;
    }
  }

}

void initialize() {

  long i,j,z;
  double r;
  
  for ( i = 0; i < MESHX; i++)
  {
    for ( j=0; j < MESHX; j++)
    {
      r= i*i + j*j;
      z= i*MESHX + j;
      if(r < 2500){
	phi_old[z] = 1.0;
      }
      else{
	phi_old[z] = 0.0;
      }
      
      mu_old[z] = Mu - deltaMu;
    }
  }
}
void boundary(double *c) {

  int i ,y ,z;

  
  
  for (i=0; i<MESHX -1; i++)
  {
    y= i*MESHX;
    z= i*MESHX + MESHX-1;
    
    //left - right
    c[y]        = c[y+1];
    c[z]  = c[z-1];
    
    // up - down
    c[i]        = c[MESHX + i];
    c[MESHX2 - MESHX + i]  = c[MESHX2 - 2*MESHX + i];
  }
  
}


void write2file ( long t) {
  int i,j,z;
  FILE *fp;
  char filename[1000];
  
  sprintf(filename,"phi_%ld.dat",t);
  fp = fopen(filename,"w");
  for ( i = 0; i < MESHX; i++)
  {
    for ( j=0; j < MESHX; j++)
    {
      
      z= i*MESHX + j;
      fprintf(fp,"%d %d %le\n",i,j, phi_new[z]);
      
    }
    fprintf(fp,"\n");
  }
  fclose(fp);
 /* 
  sprintf(filename,"chempot_%ld.dat",t);
  fp = fopen(filename,"w");
  for ( i = 0; i < MESHX; i++)
  {
    for ( j=0; j < MESHX; j++)
    {
      
      z= i*MESHX + j;
      fprintf(fp,"%d %d %le\n", i, j, mu_new[z]);
      
    }
	fprintf(fp,"\n");
  }
  fclose(fp);
  
  sprintf(filename,"conc_%ld.dat",t);
  fp = fopen(filename,"w");
  for ( i = 0; i < MESHX; i++)
  {
    for ( j=0; j < MESHX; j++)
    {
      
      z= i*MESHX + j;
      fprintf(fp,"%d %d %le\n",i,j, conc[z]);
      
    }
	fprintf(fp,"\n");
  }
  fclose(fp);  
 */
}
