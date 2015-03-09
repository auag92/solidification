#include <stdio.h>
#include <stdlib.h>

#define MESHX (128)
#define MESHX2 (MESHX*MESHX)
#define deltax (1)
#define inv_deltax (1/deltax)
#define inv_deltax2 (inv_deltax*inv_deltax)
#define deltat (.01)
#define re (100)
#define inv_re (1/re)
#define ntimesteps (100000)
#define Uw (1) // wall velocity

double P_old[MESHX2], P_now[MESHX2];
double u_old[MESHX2], u_now[MESHX2];
double v_old[MESHX2], v_now[MESHX];


double dw_dt, dxi_dx, dxi_dy, dw_dx, dw_dy, d2w_dx2, d2w_dy2;
void initialize();
void boundary();
void update();
void RHS_fn();
main(){

  int i,j,z,t=0;
  for(t=1; t<ntimesteps; t++){
    for(i=1; i<MESHX; i++){
      for(j=1; j<MESHX; j++){

        du_dt = -1*dp_dx + mu*lap_u[z] - u_old[z]*du_dx - v_old[z]*du*dy;
        dv_dt = -1*dp_dy + mu*lap_v[z] - u_old[z]*dv_dx - v_old[z]*dv*dy;
        u_new = u_old + deltat * u_old[z];
        v_new = v_old + deltat * v_old[z];
      }
    }
  }
  boundary(xi_now, w_now)
  update();

}


void RHSS_fn(){
  int i,j,x,z;

  for (i=1; i<MESHX-1; i++){
    for (j=1; j<MESHX-1; j++){
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
      fn[z] = Mu*(dlapU_dx + dlapV_dy) - advctn;
    }
  }

}
