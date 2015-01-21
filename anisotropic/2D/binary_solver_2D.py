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
Dab = 0.06

ntimesteps = 50
saveT = 1000
deltaMu = 0.3
Mu = 1.0

CONSTANT_X_0 = 0.1
CONSTANT_X_1 = 0.5

wavenumber = 4
amp = 0.1

def zero_init(n):
    return [0 for i in range(n)]

phi_old = zero_init(MESHX2)
phi_new = zero_init(MESHX2)
mu_new = zero_init(MESHX2)
mu_old = zero_init(MESHX2)
lap_phi = zero_init(MESHX2)
lap_mu = zero_init(MESHX2)
conc = zero_init(MESHX2)
dphi_now = zero_init(4*MESHX)
dphi_next = zero_init(4*MESHX)
inv_deltax2 = (1.0/deltax2)
inv_deltax = (1.0/deltax)

def update():
    for i in range(1, MESHX-1):
        for j in range(1, MESHX-1):
            z= i*MESHX + j;
            phi_old[z]=phi_new[z];
            mu_old[z]=mu_new[z];

def laplacian(f, lap):
    for i in range(1, MESHX-1):
        for j in range(1, MESHX-1):
            z = i*MESHX + j
            coeff = (f[z-1] + f[z+1] + f[z+MESHX] + f[z-MESHX]  -4.0*f[z])
            lap[z] = coeff * inv_deltax2

def concentration( ):
  for i in range(1, MESHX-1):
      for j in range(1, MESHX-1):
          z= i*MESHX + j
          p = phi_new[z]
          u = mu_new[z]
          h = p*p*(3-2*p)
          conc[z] = u*(1-h) + (h)*K*u

def initialize():
    for i in range(1, MESHX-1):
        for j in range(1, MESHX-1):
            r = i*i + j*j
            z = i*MESHX + j
            if (r < 2500.0):
                phi_old[z] = 1.0
            else:
                phi_old[z] = 0.0

            mu_old[z] = Mu - deltaMu


def boundary(c):
    for i in range(MESHX-1):
        y= i*MESHX
        z= i*MESHX + MESHX-1

        #left - right
        c[y]  = c[y+1]
        c[z]  = c[z-1]

        # up - down
        c[i]        = c[MESHX + i]
        c[MESHX2 - MESHX + i]  = c[MESHX2 - 2*MESHX + i]

def grad_phi(i, d_phi):
    if (i == MESHX -1 ):
        for j in range(1, MESHX-1):
            z = i * MESHX + j
            d_phi[2*MESHX+j] = (phi_old[z] - phi_old[z-MESHX])*inv_deltax
            d_phi[3*MESHX+j] = (phi_old[z+1] - phi_old[z-1] + phi_old[z+1-MESHX] - phi_old[z-1-MESHX])*0.25*inv_deltax
    else:
        for j in range(1, MESHX-1):
            z = i * MESHX + j
            d_phi[j] = (phi_old[z] - phi_old[z-1])*inv_deltax
            d_phi[MESHX+j] = (phi_old[z+MESHX] - phi_old[z-MESHX] + phi_old[z-1+MESHX] - phi_old[z-1-MESHX])*0.25*inv_deltax
            d_phi[2*MESHX+j] = (phi_old[z] - phi_old[z-MESHX])*inv_deltax
            d_phi[3*MESHX+j] = (phi_old[z+1] - phi_old[z-1] + phi_old[z+1-MESHX] - phi_old[z-1-MESHX])*0.25*inv_deltax
        z = i * MESHX + j
        d_phi[j] = (phi_old[z] - phi_old[z-1])*inv_deltax
        d_phi[MESHX+j] = (phi_old[z+MESHX] - phi_old[z-MESHX] + phi_old[z-1+MESHX] - phi_old[z-1-MESHX])*0.25*inv_deltax



def dqdx(phi_x, phi_y):
    ans = 0

    phi_x4 = phi_x**4
    phi_y4 = phi_y**4
    phi_x2 = phi_x**2
    phi_y2 = phi_y**2
    phi_xy2 = phi_x2 + phi_y2

    if ((phi_x2> 1e-15) and (phi_y2> 1e-15)):
        ans = 2*E*G*phi_x*(Dab*(-4*phi_x**4 - 4*phi_y**4 + 3*(phi_xy2)**2) - (phi_xy2)**2)*(Dab*(-4*phi_x**4 - 4*phi_y**4 + 3*(phi_xy2)**2) - 16*Dab*(-phi_x**4 + phi_x**2*(phi_xy2) - phi_y**4) - (phi_xy2)**2)/(phi_xy2)**4
    return ans

def div_phi(i):
    x_now = dqdx(dphi_now[i], dphi_now[i+MESHX])
    x_next = dqdx(dphi_now[i+1], dphi_now[i+1+MESHX])
    y_now = dqdx(dphi_now[i+2*MESHX], dphi_now[i+3*MESHX])
    y_next = dqdx(dphi_next[i+2*MESHX], dphi_next[i+3*MESHX])
    ans = ((x_next - x_now) + ( y_next - y_now))*inv_deltax

    return ans

def fnupdate():
    for i in range(MESHX):
        dphi_now[i] = dphi_next[i]
        dphi_now[MESHX+i] = dphi_next[MESHX+i]
        dphi_now[2*MESHX+i] = dphi_next[2*MESHX+i]
        dphi_now[3*MESHX+i] = dphi_next[3*MESHX+i]

def write2file (t):
    f = open('/home/apaar/codes/project_solidification/anisotropic/2D/datafiles2/phi'+ str(t) + '.dat','w')

    for i in range(MESHX):
        for j in range(MESHX):
            z= i*MESHX + j
            f.write('{} {} {}\n'.format(i, j , phi_new[z]))
        f.write('\n');
    f.close()


initialize()
boundary(phi_old)
boundary(mu_old)

for t in range(ntimesteps):

    laplacian(phi_old, lap_phi)
    laplacian(mu_old, lap_mu)
    grad_phi(1, dphi_now)

    for i in range(1, MESHX-1):
        grad_phi(i+1, dphi_next)
        for j in range(1, MESHX-1):
            z= i*MESHX + j
            p = phi_old[z]
            #//Gamma = 2*G*lap_phi[z]
            Gamma = div_phi(j)
            drv_frce = (mu_old[z] - Mu)*(K-1)*(mu_old[z])*6*p*(1-p)
            alln_chn = E*Gamma - (G/E)*18.0*(p)*(1.0-p)*(1.0-2.0*p)
            dp = deltat*(alln_chn + drv_frce)/(tau*E)
            phi_new[z] = p + dp
            du = deltat*M*lap_mu[z] - (K-1)*mu_old[z]*6*p*(1-p)*dp
            mu_new[z] = mu_old[z]  + du/(1+(K-1)*p*p*(3-2*p))
        fnupdate()

    boundary(phi_new)
    boundary(mu_new)
    #concentration()
    update()
    if ((t % saveT) is 0):
      write2file(t)

print("Yay! We are done.\n")
