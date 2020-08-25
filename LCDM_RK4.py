import numpy as np
from matplotlib import pyplot as plt

#Constants
pi=np.pi
pi2 = 2*np.pi
G=6.67e-11 #N m^2/kg^2
c = 9.4608e21 #m/ million years
lambd = 1e-52 #m^-1
ttoday = 13800 #million years

H0 = 6.87e-5 #1/Million years
rhocrit = 3*H0**2/8/pi/G
omega_lam = .6911
omega_m = .3089
omega_rad = 8.975e-5


# total Energy density
def rho(a,w):
    return rhocrit*(rho_cdm(a) + rho_b(a))

# photon energy density
def rho_photon(a):
    return rhophoton0 * a**4


#Lambda dark energy
def rho_lam(a,w):
    return rho_lam0 * a**(3*(1+w))

def w(a):
    return -1

#cold dark matter
def rho_cdm(a):
    return rhocdm0 * a**3

#baryon
def rho_b(a):
    return rhob0 * a**3

#curvature
def rho_k(a):
    return rhok0 * a**2

#massless neutrinos





#Friedmann equation
def adot(a,w):
    return a * H0 * (omega_lam + omega_m/a**3 + omega_rad/a**4)**(1/2)

#Runge Kutta
def rk4(t0, tf, a0, h):   #initial time,final time, initial condition, step
    n = (int)((tf)/h)
    t=np.linspace(0, tf, n)
    a = np.zeros([n])
    j = len(t)-1
    a[j] = a0
    while a[j]>0 and t[j]>0 and j>0:
        k1 = h * adot(a[j],0)
        k2 = h * adot(a[j]+k1/2,0)
        k3 = h * adot(a[j]+k2/2, 0)
        k4 = h * adot(a[j]+k3, 0)
        a[j-1] = a[j] - (k1/6 + k2/3 + k3/3 + k4/6)
        j -= 1
    return a
    return t





plt.plot(rk4(0,ttoday,1,1e-1))
plt.show()