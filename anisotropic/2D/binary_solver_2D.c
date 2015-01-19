#include "stdio.h"
#include "math.h"
#include <sys/stat.h>

#define MESHX 200
#define MESHX2 (MESHX*MESHX)
#define deltat (0.02)
#define deltax (1)
#define deltax2 (deltax*deltax)

#define K (0.5)             /*Partition Coefficient*/
#define G (1.0)             /*Surface Energy*/
#define M (1.0)             /*Mobility*/
#define E (4.0)             /*Relaxation factor - dimensions of length [m]*/
#define tau (1.0)
#define Dab (1)

#define ntimesteps (1)
#define saveT (1000)
#define deltaMu (0.3)
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
double dphi_old[MESHX], dphi_new[MESHX];
double inv_deltax2 = (1.0/deltax2);
double inv_deltax = (1.0/deltax);


void initialize();
void boundary(double *c);
void write2file ( long t);
void update();
void laplacian(double *f, double *lap);
void concentration();
double aniso_x( long i, long j);
double aniso_y( long i, long j);
void fnupdade();

void main() {
  long i, j, z, t;
  double p,dp,du;

  double Gamma;
  initialize();
  boundary(phi_old);
  boundary(mu_old);

  for (t=0; t < ntimesteps; t++) {
    //laplacian(phi_old,     lap_phi);
    laplacian(mu_old,     lap_mu);


    for (i=1; i < (MESHX-1); i++) {

      for (j=1; j < (MESHX-1); j++){

        z= i*MESHX + j;
        p = phi_old[z];
        Gamma =
        //Gamma = 2*G*lap_phi[z];
        dp = deltat*(E*Gamma - (G/E)*18.0*(p)*(1.0-p)*(1.0-2.0*p) + (mu_old[z] - Mu)*(K-1)*(mu_old[z])*6*p*(1-p))/(tau*E);
        phi_new[z] = p + dp;
        du = deltat*M*lap_mu[z] - (K-1)*mu_old[z]*6*p*(1-p)*dp;
        mu_new[z] = mu_old[z]  + du/(1+(K-1)*p*p*(3-2*p));

        x_now=x_next;
      }
      fnupdade();
    }
    boundary(phi_new);
    boundary(mu_new);
    //concentration();
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
      lap[z] = (f[z-1] + f[z+1] + f[z+MESHX] + f[z-MESHX]  -4.0*f[z])*inv_deltax2;
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
      if(r < 2500.0){
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
void grad_phi(double *phi, double *d_phi){

 	long i,z;

 	for(i=1; i<MESHX-1; i++){
 		z = MESHX + i;
 		d_phi[i] = (phi_old[z] - phi_old[z-1])*inv_deltax;
        	d_phi[MESHX+i] = (phi_old[z+MESHX] - phi_old[z-MESHX] + phi_old[z-1+MESHX] - phi_old[z-1-MESHX])*0.25*inv_deltax;
        	d_phi[2*MESHX+i] = (phi_old[z] - phi_old[z-MESHX])*inv_deltax;
 		d_phi[3*MESHX+i] = (phi_old[z+1] - phi_old[z-1] + phi_old[z+1-MESHX] - phi_old[z-1-MESHX])*0.25*inv_deltax;
 	}
 	for(i=1; i<MESHX-1; i++){

 		d_phi[i] = (phi_old[z] - phi_old[z-1])*inv_deltax;
        	d_phi[MESHX+i] = (phi_old[z+MESHX] - phi_old[z-MESHX] + phi_old[z-1+MESHX] - phi_old[z-1-MESHX])*0.25*inv_deltax;
        	d_phi[2*MESHX+i] = (phi_old[z] - phi_old[z-MESHX])*inv_deltax;
 		d_phi[3*MESHX+i] = (phi_old[z+1] - phi_old[z-1] + phi_old[z+1-MESHX] - phi_old[z-1-MESHX])*0.25*inv_deltax;
 	}


}
double dqdx( double phi_x, double phi_y){

	double a, ans, phi_x2, phi_x4, phi_y2, phi_y4, inv_phi2;
	double ans;
	long z;


        phi_x4 = phi_x *phi_x *phi_x *phi_x;
        phi_y4 = phi_y *phi_y *phi_y *phi_y;
        phi_x2 = phi_x *phi_x;
        phi_y2 = phi_y *phi_y;

	inv_phi = 1/(phi_x2+phi_y2);

        a = G*(1 - Dab*(4*(phi_x4 + phi_y4) - 3)*inv_phi*inv_phi);
        ans= 2 * a * phi_x * G * (1 - Dab*(16.0*phi_x2 + (-12.0*(phi_x4+phi_y4)+9.0)*inv_phi)*inv_phi);
        return ans;
}

double div_phi(){
	double ans;
	double x_next, x_now;
	double y_next, y_now;
	x_next = 0;
	ans = ((x_next - x_now) + ( y_next[j]-y_now[j]))*inv_deltax;
	return ans;
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
}
