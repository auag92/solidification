#include "stdio.h"
#include "math.h"
#include <sys/stat.h>

#define MESHX 200
#define deltat (0.02)
#define deltax (1.0)
#define deltax2 (deltax*deltax)


#define K (0.5)             /*Partition Coefficient*/
#define G (1.0)             /*Surface Energy*/
#define M (1.0)             /*Mobility*/
#define E (4.0)            /*epsilon - dimensions of length [m]*/
#define tau (1.0)

#define ntimesteps (500000)
#define saveT (1000)

#define deltaMu (0.1)
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

double phi_new[MESHX], phi_old[MESHX];
double mu_new[MESHX], mu_old[MESHX];
double lap_phi[MESHX], lap_mu[MESHX];
double conc[MESHX];

double inv_deltax2 = (1.0/deltax2);

void initialize(double *c);
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
double compute_surface_energy (double *c);

void main() {
  long i, t;
  double p,dp,du;
  double surface_energy;
  
  initialize(phi_old);
  boundary(phi_old);
  boundary(mu_old);

  for (t=0; t < ntimesteps; t++) {

    laplacian(phi_old,     lap_phi);
    laplacian(mu_old,     lap_mu);

    for (i=0; i < (MESHX); i++) {

      p = phi_old[i];
      dp = deltat*(G*E*lap_phi[i] - (G/E)*18.0*(p)*(1.0-p)*(1.0-2.0*p) + (mu_old[i] - Mu)*(K-1)*(mu_old[i])*6*p*(1-p))/(tau*E);
      phi_new[i] = p + dp;
      du = deltat*M*lap_mu[i] - (K-1)*mu_old[i]*6*p*(1-p)*dp;
      mu_new[i] = mu_old[i]  + du/(1+(K-1)*p*p*(3-2*p));

    }

    boundary(phi_new);
    boundary(mu_new);
    concentration();
    update();
    if((t%saveT) ==0) {
      write2file(t);
      
      surface_energy = compute_surface_energy(phi_old);
      printf("surface_energy=%le\n", surface_energy);
    }
      

  }
  printf("Yay! We are done.\n");
}

void update() {
  long i;
  for (i=0; i < MESHX; i++) {
    phi_old[i]=phi_new[i];
    mu_old[i]=mu_new[i];

  }
}

void laplacian(double *f, double *lap) {
  long i;
  for (i=1; i < (MESHX-1); i++) {
    lap[i] = (f[i-1] + f[i+1] -2.0*f[i])*inv_deltax2;
  }
}

void concentration(){

  double p,u,h;
  int i;
  for ( i = 0; i < MESHX; i++)
  {
    p = phi_new[i];
    u = mu_new[i];
    h = p*p*(3-2*p);
    conc[i] = u*(1-h) + (h)*K*u;
  }

}

void initialize(double *c) {

  long i;
  #ifdef STEP_PROFILE
  for (i=0; i < MESHX/4; i++) {
    c[i] = 1.0;
  }
  for (i=MESHX/4; i < MESHX; i++) {
    c[i] = 0.0;
  }
  #endif

  #ifdef SINE_PROFILE


  for (i=1; i < MESHX-1; i++) {
    c[i]= 0.5 + amp*sin(2.0*M_PI*(i-0.5)*wavenumber/(MESHX-1));
  }
  #endif

  for (i=1; i < MESHX-1; i++) {
    mu_old[i] = Mu - deltaMu;
  }
}
void boundary(double *c) {

  #ifdef DIRICHLET
  c[0]        = CONSTANT_X_0;
  c[MESHX-1]  = CONSTANT_X_1;

  #endif

  #ifdef NEUMANN
  c[0]        = c[1];
  c[MESHX-1]  = c[MESHX-2];

  #endif

  #ifdef PERIODIC
  c[0]            = c[MESHX-2];
  c[MESHX-1]      = c[1];

  #endif
}

double compute_surface_energy (double *c) {
  long i;
  double integral=0.0;
  
  for (i=2; i < MESHX-2; i++) {
    integral += 0.5*(c[i]*c[i]*(1.0-c[i])*(1.0-c[i]) + c[i+1]*c[i+1]*(1.0-c[i+1])*(1.0-c[i+1])); 
  }
  return integral*deltax;
}



void write2file ( long t) {
  long i;
  FILE *fp;
  char filename[1000];

  sprintf(filename,"phi_%ld.dat",t);
  fp = fopen(filename,"w");
  for (i=0; i < (MESHX); i++) {
    fprintf(fp,"%le %le\n", i*deltax, phi_new[i]);
  }
  fclose(fp);
  sprintf(filename,"ChemPot_%ld.dat",t);
  fp = fopen(filename,"w");
  for (i=0; i < (MESHX); i++) {
    fprintf(fp,"%le %le\n", i*deltax, mu_new[i]);
  }
  fclose(fp);

  sprintf(filename,"Concentration_%ld.dat",t);
  fp = fopen(filename,"w");
  for (i=0; i < (MESHX); i++) {
    fprintf(fp,"%le %le\n", i*deltax, conc[i]);
  }
  fclose(fp);

}
