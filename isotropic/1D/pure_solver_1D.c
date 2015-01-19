#include"stdio.h"
#include"math.h"

#define SOLIDIFICATION

#define MESHX 500
#define deltat (0.03)
#define deltax (1.0)
#define deltax2 (deltax*deltax)
#define D (1.0)
#define K (1.0)             /*Condutivity*/
#define L (1.0)             /*Latent heat*/
#define Cv (1.0)            /*Specific heat*/
#define M (1.0)

#define sigma (1.0)  /*Surface energy*/
#define W (10.0)     /*Interface-width*/
#define H (3.0*sigma*4.4/W)  /*Height of potential*/
#define kappa (3.0*sigma*W/4.4)


#define ntimesteps (5000000)
#define saveT (10000)
#define DELTAT 0.1
#define Tm (1.0)

// // #define DIRICHLET
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


double c_new[MESHX], c_old[MESHX];
double T_new[MESHX], T_old[MESHX];
double mu[MESHX],    lap_mu[MESHX], lap_c[MESHX], lap_T[MESHX], lap4_c[MESHX], g[MESHX];
double inv_deltax2 = (1.0/deltax2);

void initialize(double *c);
void boundary(double *c);
void write2file (double *c, long t);
void update();
void compute_chemical_potential(double *c, double *mu);
void laplacian(double *f, double *lap);
void free_energy(double *c, double *g);
double compute_surface_energy(double *c);

void main() {
  long i, t;

  initialize(c_old);
  boundary(c_old);
  boundary(T_old);

  for (t=0; t < ntimesteps; t++) {


   laplacian(c_old,     lap_c);
   laplacian(T_old,     lap_T);
   compute_chemical_potential(c_old, mu); //Computes (df_o/deta)

   for (i=0; i < (MESHX); i++) {
      c_new[i] = c_old[i] + (deltat)*M*(-H*mu[i] + 2.0*kappa*lap_c[i] -6.0*L*(T_old[i] - Tm)*c_old[i]*(1.0-c_old[i])/Tm);

      T_new[i] = T_old[i] + (deltat/Cv)*(K*lap_T[i])
               + (L/Cv)*6.0*c_old[i]*(1.0-c_old[i])*(c_new[i]-c_old[i]);

    }
#endif
    boundary(c_new);
    boundary(T_new);
    update();
    if((t%saveT) ==0) {
      write2file(c_new,t);
    }
  }
  double surface_energy;
  printf("Yay! We are done.\n");
}

void update() {
  long i;
  for (i=0; i < MESHX; i++) {
    c_old[i]=c_new[i];
#ifdef SOLIDIFICATION
    T_old[i]=T_new[i];
#endif
  }
}
void compute_chemical_potential(double *c, double *mu) {
  long i;
  for (i=0; i < MESHX; i++) {
    mu[i] = 2.0*c[i]*(1.0 -3.0*c[i] + 2.0*c[i]*c[i]);
  }
}
void free_energy(double *c, double *g) {
  long i;

#ifdef SOLIDIFICATION
  for (i=0; i < MESHX; i++) {
    g[i] = H*c[i]*c[i]*(1.0-c[i])*(1.0-c[i]);
  }
#endif
}

void laplacian(double *f, double *lap) {
  long i;
  for (i=1; i < (MESHX-1); i++) {
    lap[i] = (f[i-1] + f[i+1] -2.0*f[i])*inv_deltax2;
  }
//   lap[0]       = 0.0;
//   lap[MESHX-1] = 0.0;
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

#ifdef SOLIDIFICATION
  for (i=1; i < MESHX-1; i++) {
    c[i]= 0.5 + amp*sin(2.0*M_PI*(i-0.5)*wavenumber/(MESHX-1));
  }
#endif
#endif

#ifdef SOLIDIFICATION
  for (i=1; i < MESHX-1; i++) {
    T_old[i] = Tm - DELTAT;
  }
#endif
}
void boundary(double *c) {

#ifdef SOLIDIFICATION
  #ifdef DIRICHLET
  c[0]        = CONSTANT_X_0;
  c[MESHX-1]  = CONSTANT_X_1;

  lap_c[0]        = 0.0;
  lap_c[MESHX-1]  = 0.0;
  #endif

  #ifdef NEUMANN
  c[0]        = c[1];
  c[MESHX-1]  = c[MESHX-2];

  lap_c[0]        = 0.0;
  lap_c[MESHX-1]  = 0.0;
  #endif

  #ifdef PERIODIC
  c[0]            = c[MESHX-2];
  c[MESHX-1]      = c[1];

  lap_c[0]        = lap_c[MESHX-2];
  lap_c[MESHX-1]  = lap_c[1];
  #endif
#endif
}

void write2file (double *c, long t) {
  long i;
  FILE *fp;
  char filename[1000];

#ifdef SOLIDIFICATION
  sprintf(filename,"phi_%ld.dat",t);
  fp = fopen(filename,"w");
  for (i=0; i < (MESHX); i++) {
    fprintf(fp,"%le %le\n", i*deltax, c[i]);
  }
  fclose(fp);
  sprintf(filename,"Temperature_%ld.dat",t);
  fp = fopen(filename,"w");
  for (i=0; i < (MESHX); i++) {
    fprintf(fp,"%le %le\n", i*deltax, T_new[i]);
  }
  fclose(fp);
#endif
}
