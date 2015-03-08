
for (i=1; i<MESHX-1; i++)
{
  for (j=1; j<MESHX-1; j++){
    z = i*MESHX + j;
    dxi_dx = xi_old[z+1] - xi_old[z-1];
    dxi_dy = xi_old[z+MESHX] - xi_old[z-MESHX];
    dw_dx = w_old[z+1] - w_old[z-1];
    dw_dy = w_old[z+MESHX] - w_old[z-MESHX];
    lap_w = w_old[z+MESHX]+w_old[z-MESHX]+w_old[z+1]+w_old[z-1]-4*w_old[z];
    dw_dt = inv_deltax2*(0.25*(-1*dxi_dy*dw_dx + dxi_dx*dw_dy) + inv_re*lap_w);
    w_now[z] = w_old[z] + deltat * dw_dt;



  }
