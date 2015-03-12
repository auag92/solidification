#include <stdio.h>
#include <stdlib.h>
#include "mg_solver.c"
#include "math.h"
#define pmesh 65
#define pmesh2 (pmesh*pmesh)
#define MESHX (pmesh + 1)
#define MESHX2 (MESHX*MESHX)
#define deltax (1)
#define inv_deltax (1/deltax)
#define inv_deltax2 (inv_deltax*inv_deltax)
#define deltat (.01)
#define ntimesteps (1)
#define saveT (100)

#define Mu (1)
#define re (100)
#define inv_re (1/re)
#define Uw (1) // wall velocity

double P[2*pmesh2];
double u_old[MESHX2], u_now[MESHX2];
double v_old[MESHX2], v_now[MESHX];
double rhs_fn[2*pmesh2];

void initialize();
void boundary1();
void update1(double *old, double *new,int M);
void laplacian(double *f, double *lap, int M);
void RHS_fn(double *u, double *v, double *lap_u, double *lap_v, int M);
void printArray(double *c, int M);
void write2file (long t, double *u, double *v, int M);

main(){

  int i,j,z,x,t=0;
  double du_dx,du_dy,dv_dx,dv_dy,dp_dx,dp_dy,du_dt,dv_dt;
  double lap_u[MESHX2], lap_v[MESHX2];
  initialize();
  boundary1();
  for(t=0; t<ntimesteps; t++){
    laplacian(u_old,lap_u, MESHX);
    laplacian(v_old,lap_v, MESHX);
    RHS_fn(u_old,v_old,lap_u,lap_v,MESHX);
    multigrid(P, rhs_fn);
    printArray(P, pmesh);
    for(i=1; i<MESHX-1; i++){
      for(j=1; j<MESHX-1; j++){
	
	z = i*MESHX +j;
	x = (i-1)*(pmesh)+ (j-1); 
	du_dx = 0.5*inv_deltax*(u_old[z+1] - u_old[z-1]);
	dv_dx = 0.5*inv_deltax*(v_old[z+MESHX] - v_old[z-MESHX]);
	dp_dx = 0.5*inv_deltax*(P[x+1]+P[x+pmesh+1]-P[x+pmesh]-P[x]);
	dp_dy = 0.5*inv_deltax*(P[x+pmesh]+P[x+pmesh+1]-P[x+1]-P[x]);
        du_dt = Mu*lap_u[z] - dp_dx - (u_old[z]*du_dx + v_old[z]*du_dy);
        dv_dt = Mu*lap_v[z] - dp_dy - (u_old[z]*dv_dx + v_old[z]*dv_dy);
        u_now[z] = u_old[z] + deltat * u_old[z];
        v_now[z] = v_old[z] + deltat * v_old[z];
      }
    }
    boundary1();
    update1(v_old, v_now, MESHX);
    update1(u_old, u_now, MESHX);
    if((t%saveT) ==0) {
      write2file(t,u_old,v_old, MESHX);
    }    
  }
  printf("Yay! We are done.\n");
}
void update1(double *old, double *new,int M) {
  long i, j, z;
  for (i=0; i < M; i++) {
    for (j=0; j < M; j++){
      z= i*M + j;
      old[z]=new[z];
      old[z]=new[z];      
    }
  }
}
  
void laplacian(double *f, double *lap, int M) {
  long i,j,z;
    
  for (i=1; i< M -1; i++)
  {
    for (j=1; j< M -1; j++)
    {
      
      z= i*M + j;
      lap[z] = (f[z-1] + f[z+1] + f[z+M] + f[z-M] -4.0*f[z])*inv_deltax2;

    }
  }
}
void initialize() {
  long i,j,z;
    
  for (i=0; i< MESHX-1; i++)
  {
    for (j=0; j< MESHX -1; j++)
    {      
      v_old[z] = 0.0;
      u_old[z] = 0.0;
    }
  }
  
  for (i=0; i< pmesh-1; i++)
  {
    for (j=0; j< pmesh -1; j++)
    {      
      z= i*pmesh + j;
      P[z] = 0.0;
      rhs_fn[z] = 0.0;
    }
  }
}
void boundary1() {

  int i ,y ,z; 
  for (i=0; i<MESHX; i++)
  {
    u_now[i] = Uw;
    v_now[i] = 0.0;
    u_now[i*MESHX] = 0.0;
    v_now[i*MESHX] = 0.0;
    u_now[i*MESHX+MESHX-1] = 0.0;
    v_now[i*MESHX+MESHX-1] = 0.0;
    u_now[MESHX2-MESHX+i] = 0.0;
    v_now[MESHX2-MESHX+i] = 0.0;
  }  
}

void RHS_fn(double *u, double *v, double *lap_u, double *lap_v, int M){
  int i,j,x,z;
  double du_dx,du_dy,u_avg,d2u_dx2,d2u_dy2,d2u_dxdy,dv_dx,dv_dy,v_avg,d2v_dy2,d2v_dydx;
  double advctn, dlapU_dx,dlapV_dy;
 
  for (i=1; i<M-2; i++){
    for (j=1; j<M-2; j++){
      
      z = i*M + j;
      x = i*(M-1) + j;
      dlapU_dx =0.5*inv_deltax*(lap_u[z+1]-lap_u[z]+lap_u[z+M+1]-lap_u[z+M]);
      dlapV_dy =0.5*inv_deltax*(lap_v[z+M]-lap_v[z]+lap_v[z+1+M]-lap_v[z+1]);

      du_dx = 0.5*inv_deltax*(u[z+1]-u[z]+u[z+M+1]-u[z+M]);
      du_dy = 0.5*inv_deltax*(u[z+M]-u[z]+u[z+1+M]-u[z+1]);
      u_avg = 0.25*(u[z]+u[z+1]+u[z+M]+u[z+M+1]);
      d2u_dx2 = 0.5*inv_deltax2*(u[z-1]-u[z]-u[z+1]+u[z+2]+u[z-1+M]-u[z+M]-u[z+1+M]+u[z+2+M]);
      d2u_dxdy = inv_deltax2*(u[z+M+1]-u[z+1]-u[z+M]+u[z]);

      dv_dx = 0.5*inv_deltax*(v[z+1]-v[z]+v[z+M+1]-v[z+M]);
      dv_dy = 0.5*inv_deltax*(v[z+M]-v[z]+v[z+1+M]-v[z+1]);
      v_avg = 0.25*(v[z]+v[z+1]+v[z+M]+v[z+M+1]);
      d2v_dy2 = 0.5*inv_deltax2*(v[z-M]-v[z]-v[z+M]+v[z+2*M]+v[z-M+1]-v[z+1]-v[z+M+1]+v[z+2*M+1]);
      d2v_dydx = inv_deltax2*(v[z+M+1]-v[z+M]-v[z+1]+v[z]);

      advctn = du_dx*du_dx+dv_dy*dv_dy+2*(dv_dx*du_dy)+v_avg*(d2v_dy2+d2u_dxdy)+u_avg*(d2u_dx2+d2v_dydx);
      rhs_fn[x] = Mu*(dlapU_dx + dlapV_dy) - advctn;
    }
  }
  printArray(rhs_fn, pmesh);

}
void printArray(double *c, int M){
  long i,j,z;
    
  for (i=0; i< M-1; i++)
  {
    for (j=0; j< M -1; j++)
    {  
      z= i*M+j;
      printf("%le ",c[z]);
    }
    printf("\n");
  }
  
}
void write2file (long t, double *u, double *v, int M) {
  int i,j,z;
  FILE *fp;
  char filename[1000];
  sprintf(filename,"./datafiles/velocity_%ld.dat",t);
  fp = fopen(filename,"w");
  if (fp) {
    for ( i = 0; i < M; i++)
    {
      for ( j=0; j < M; j++)
      {
        z= i*M + j;
        fprintf(fp,"%d %d %le %le\n",i,j, u[z],v[z]);
      }
      fprintf(fp,"\n");
    }
    fclose(fp);
  }
  else {
    printf("#Error:%s could not be opened to write",filename);
    exit(1);
    }
  }