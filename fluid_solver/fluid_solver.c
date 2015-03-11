#include <stdio.h>
#include <stdlib.h>

#define pmesh 128
#define MESHX (pmesh + 3)
#define MESHX2 (MESHX*MESHX)
#define deltax (1)
#define inv_deltax (1/deltax)
#define inv_deltax2 (inv_deltax*inv_deltax)
#define deltat (.01)
#define ntimesteps (100000)

#define Mu (1)
#define re (100)
#define inv_re (1/re)
#define Uw (1) // wall velocity

double P[MESHX2], Pstr[pmesh];
double u_old[MESHX2], u_now[MESHX2];
double v_old[MESHX2], v_now[MESHX];
double rhs_fn[pmesh];

void initialize();
void boundary();
void update();
void laplacian(double *f, double *lap);
void RHS_fn(double *u, double *v, double *lap_u, double *lap_v, int M);

main(){

  int i,j,z,t=0;
  double du_dx,du_dy,dv_dx,dv_dy,dp_dx,dp_dy,du_dt,dv_dt;
  double lap_u[MESHX2], lap_v[MESHX2];
  initialize();
  for(t=0; t<ntimesteps; t++){
    laplacian(*u_old,*lap_u, MESHX2);
    laplacian(*v_old,*lap_v, MESHX2);
    RHS_fn(*u_old,*v_old,*lap_u,*lap_v,*p_old,MESHX2);
    multigrid(*Pstr, *rhs_fn, pmesh);
    for(i=1; i<MESHX-1; i++){
      for(j=1; j<MESHX-1; j++){
	du_dx = 0.5*inv_deltax*(u_old[z+i] - u_old[z-i]);
	dv_dx = 0.5*inv_deltax*(v_old[z+MESHX] - v_old[z-MESHX]);
	dp_dx = inv_deltax*(P[z+i]-P[z-i]);
	dp_dy = inv_deltax*(P[z+MESHX]-P[z-MESHX]); 
        du_dt = mu*lap_u[z] - dp_dx - (u_old[z]*du_dx + v_old[z]*du_dy);
        dv_dt = mu*lap_v[z] - dp_dy - (u_old[z]*dv_dx + v_old[z]*dv_dy);
        u_now[z] = u_old[z] + deltat * u_old[z];
        v_now[z] = v_old[z] + deltat * v_old[z];
      }
    }
  }
  boundary(v_now, u_now);
  update(v_old, v_now);
  update(u_old, u_now);
}
void update(double *old, double *new) {
  long i, j, z;
  for (i=0; i < MESHX; i++) {
    for (j=0; j < MESHX; j++){
      z= i*MESHX + j;
      old[z]=new[z];
      old[z]=new[z];      
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
      lap[z] = (f[z-1] + f[z+1] + f[z+MESHX] + f[z-MESHX] -4.0*f[z])*inv_deltax2;

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
// to be added
    }
  }
}
void boundary(double *c) {

  int i ,y ,z;

  
  
  for (i=0; i<MESHX -1; i++)
  {
// to be added
  }
  
}

void RHS_fn(double *u, double *v, double *lap_u, double *lap_v, int M){
  int i,j,x,z;
  double du_dx,du_dy,u_avg,d2u_dx2,d2u_dy2,d2u_dxdy,dv_dx,dv_dy,v_avg,d2v_dy2,d2v_dydx;
  double advctn;
 
  for (i=1; i<MESHX-2; i++){
    for (j=1; j<MESHX-2; j++){
      
      z = i*MESHX + j;
      x = z - MESHX -1;
      dlapU_dx =0.5*inv_deltax*(lap_u[z+i]-lap_u[z]+lap_u[z+MESHX+i]-lap_u[z+MESHX]);
      dlapV_dy =0.5*inv_deltax*(lap_v[z+meshx]-lap_v[z]+lap_v[z+i+meshx]-lap_v[z+i]);

      du_dx = 0.5*inv_deltax*(u[z+i]-u[z]+u[z+MESHX+i]-u[z+MESHX]);
      du_dy = 0.5*inv_deltax*(u[z+meshx]-u[z]+u[z+i+meshx]-u[z+i]);
      u_avg = 0.25*(u[z]+u[z+i]+u[z+MESHX]+u[z+MESHX+i]);
      d2u_dx2 = 0.25*inv_deltax2*(u[z-1]-u[z]-u[z+1]+u[z+2]+u[z-1+MESHX]-u[z+MESHX]-u[z+1+MESHX]+u[z+2+MESHX]);
      d2u_dxdy = inv_deltax2*(u[z+MESHX+i]-u[z+i]-u[z+MESHX]+u[z]);

      dv_dx = 0.5*inv_deltax*(v[z+i]-v[z]+v[z+MESHX+i]-v[z+MESHX]);
      dv_dy = inv_deltax*(v[z+meshx]-v[z]+v[z+i+meshx]-v[z+i]);
      v_avg = 0.25*(v[z]+v[z+i]+v[z+MESHX]+v[z+MESHX+i]);
      d2v_dy2 = 0.25*inv_deltax2*(v[z-MESHX]-v[z]-v[z+MESHX]+v[z+2*MESHX]+v[z-MESHX+i]-v[z+i]-v[z+MESHX+i]+v[z+2*MESHX+i]);
      d2v_dydx = inv_deltax2*(v[z+MESHX+i]-v[z+MESHX]-v[z+i]+v[z]);

      advctn = du_dx*du_dx+dv_dy*dv_dy+2*(dv_dx*du_dy)+v_avg*(d2v_dy2+d2u_dxdy)+u_avg*(d2u_dx2+d2v_dydx);
      rhs_fn[x] = Mu*(dlapU_dx + dlapV_dy) - advctn;
    }
  }

}
void write2file (long t, *phi) {
  int i,j,z;
  FILE *fp;
  char filename[1000];
  sprintf(filename,"./datafiles2/phi_%ld.dat",t);
  fp = fopen(filename,"w");
  if (fp) {
    for ( i = 0; i < MESHX; i++)
    {
      for ( j=0; j < MESHX; j++)
      {
        z= i*MESHX + j;
        fprintf(fp,"%d %d %le\n",i,j, phi[z]);
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