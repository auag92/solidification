#!/usr/bin/env python3

MESHX = 200
MESHX2 = MESHX*MESHX
deltat = 0.02
deltax = 1
deltax2  = deltax*deltax

K = 0.5         #    /*Partition Coefficient*/
G = 1.0         #    /*Surface Energy*/
M = 1.0         #    /*Mobility*/
E = 4.0         #    /*Relaxation factor - dimensions of length [m]*/
tau = 1.0
Dab = 0.0000

ntimesteps = 100000
saveT = 1000
deltaMu = 0.3
Mu = 1.0

#DIRICHLET
CONSTANT_X_0 = 0.1
CONSTANT_X_1 = 0.5
#NEUMANN

#PERIODIC
#SINE_PROFILE
#ifdef SINE_PROFILE
wavenumber = 4
amp = 0.1
#endif
#STEP_PROFILE

def laplacian(f, lap):
  for i in range(MESHX):
      for j in range(MESHX):
          z = i*MESHX + j
          coeff = (f[z-1] + f[z+1] + f[z+MESHX] + f[z-MESHX]  -4.0*f[z])
          lap[z] = coeff * inv_deltax2

def concentration( ):
  for i in range(MESHX):
      for j in range(MESHX):
          z= i*MESHX + j
          p = phi_new[z]
          u = mu_new[z]
          h = p*p*(3-2*p)
          conc[z] = u*(1-h) + (h)*K*u

def initialize():
    for i in range(MESHX):
        for j in range(MESHX):
            r = i*i + j*j;
            z = i*MESHX + j;
            if (r < 2500.0):
                phi_old[z] = 1.0;
            else:
                phi_old[z] = 0.0;

            mu_old[z] = Mu - deltaMu;


def boundary(c):
    for i in range(MESHX):
        y= i*MESHX
        z= i*MESHX + MESHX-1

        #//left - right
        c[y]  = c[y+1]
        c[z]  = c[z-1]

        #// up - down
        c[i]        = c[MESHX + i]
        c[MESHX2 - MESHX + i]  = c[MESHX2 - 2*MESHX + i]

void grad_phi(long i, double *d_phi){

 	long j,z;

  if (i == MESHX -1 ){
    for(j=1; j<MESHX-1; j++){
     z = i * MESHX + j;
     d_phi[2*MESHX+j] = (phi_old[z] - phi_old[z-MESHX])*inv_deltax;
     d_phi[3*MESHX+j] = (phi_old[z+1] - phi_old[z-1] + phi_old[z+1-MESHX] - phi_old[z-1-MESHX])*0.25*inv_deltax;
    }
  }
  else {
    for(j=1; j<MESHX-1; j++){
       z = i * MESHX + j;
       d_phi[j] = (phi_old[z] - phi_old[z-1])*inv_deltax;
       d_phi[MESHX+j] = (phi_old[z+MESHX] - phi_old[z-MESHX] + phi_old[z-1+MESHX] - phi_old[z-1-MESHX])*0.25*inv_deltax;
       d_phi[2*MESHX+j] = (phi_old[z] - phi_old[z-MESHX])*inv_deltax;
       d_phi[3*MESHX+j] = (phi_old[z+1] - phi_old[z-1] + phi_old[z+1-MESHX] - phi_old[z-1-MESHX])*0.25*inv_deltax;
     }
     z = i * MESHX + j;
     d_phi[j] = (phi_old[z] - phi_old[z-1])*inv_deltax;
     d_phi[MESHX+j] = (phi_old[z+MESHX] - phi_old[z-MESHX] + phi_old[z-1+MESHX] - phi_old[z-1-MESHX])*0.25*inv_deltax;
  }
}
