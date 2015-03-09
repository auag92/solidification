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
  dlapU_dx = inv_deltax*(lap_u[z+i] - lap_u[z]);
  dlapV_dy = inv_deltax*(lap_v[z+meshx] - lap_v[z]);
  du_dx = inv_deltax*(u[z+i] - u[z]);
  ud2u_dx2 = ????????
  du_dy = inv_deltax*(u[z+meshx] - u[z]);
  dv_dx = inv_deltax*(v[z+i] - v[z]);
  vd2v_dy2 = ????????
  dv_dy = inv_deltax*(v[z+meshx] - v[z]);
  advctn = du_dx*du_dx + dv_dy*dv_dy + 2*dv_dx*du_dy + vd2v_dy2 + ud2u_dx2;
  fn[z] = Mu*(dlapU_dx + dlapV_dy) - advctn;
}
