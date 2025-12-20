'''
Fourier stability anaysis of WENO for linear advection
   u_t + a u_x = 0
We set a = 1.

See:
    Wang, Spiteri: LINEAR INSTABILITY OF THE FIFTH-ORDER WENO METHOD
    SIAM J. Num. Anal., Vol. 45, No. 5, 2007
'''
import numpy as np
import matplotlib.pyplot as plt
import argparse

# Coefficients of 3-stage SSPRK
ark = [0.0, 3.0/4.0, 1.0/3.0]
a = 1.0

# 5th order weno reconstruction of Jiang-Shu. Gives left state at
# interface b/w u0 and up1
def weno5(um2,um1,u0,up1,up2):
    eps = 1.0e-6; gamma1=1.0/10.0; gamma2=3.0/5.0; gamma3=3.0/10.0;
    beta1 = (13.0/12.0)*(um2 - 2.0*um1 + u0)**2 + (1.0/4.0)*(um2 - 4.0*um1 + 3.0*u0)**2
    beta2 = (13.0/12.0)*(um1 - 2.0*u0 + up1)**2 + (1.0/4.0)*(um1 - up1)**2
    beta3 = (13.0/12.0)*(u0 - 2.0*up1 + up2)**2 + (1.0/4.0)*(3.0*u0 - 4.0*up1 + up2)**2

    w1 = gamma1 / (eps+beta1)**2
    w2 = gamma2 / (eps+beta2)**2
    w3 = gamma3 / (eps+beta3)**2

    u1 = (1.0/3.0)*um2 - (7.0/6.0)*um1 + (11.0/6.0)*u0
    u2 = -(1.0/6.0)*um1 + (5.0/6.0)*u0 + (1.0/3.0)*up1
    u3 = (1.0/3.0)*u0 + (5.0/6.0)*up1 - (1.0/6.0)*up2

    return (w1 * u1 + w2 * u2 + w3 * u3)/(w1 + w2 + w3)

def rhs(nu,theta):
    # Returns exp(i * l * theta)
    u = lambda l: np.exp(1j * l * theta)

    j = 0
    ujm3 = u(j-3)
    ujm2 = u(j-2)
    ujm1 = u(j-1)
    uj   = u(j)
    ujp1 = u(j+1)
    ujp2 = u(j+2)
    ujp3 = u(j+3) # Not needed for this pde

    fjmh = weno5(ujm3, ujm2, ujm1, uj, ujp1)
    fjph = weno5(ujm2, ujm1, uj, ujp1, ujp2)
    z = - nu * (fjph - fjmh) / uj
    return z

# RK1 Amplification factor
def gRK1(nu, theta):
    z = rhs(nu,theta)
    return 1 + z

# RK2 Amplification factor
def gRK2(nu, theta):
    z = rhs(nu,theta)
    return 1 + z + z**2/2

# RK3 Amplification factor
def gRK3(nu, theta):
    z = rhs(nu,theta)
    return 1 + z + z**2/2 + z**3/6

# RK4 Amplification factor
def gRK4(nu, theta):
    z = rhs(nu,theta)
    return 1 + z + z**2/2 + z**3/6 + z**4/24

ntheta = 1000
theta = np.linspace(0, np.pi, ntheta)
gvalues1 = np.zeros(ntheta, dtype=complex)
gvalues2 = np.zeros(ntheta, dtype=complex)
gvalues3 = np.zeros(ntheta, dtype=complex)
gvalues4 = np.zeros(ntheta, dtype=complex)

# CFL
nu = 0.1

for k in range(ntheta):
    gvalues1[k] = gRK1(nu, theta[k])
    gvalues2[k] = gRK2(nu, theta[k])
    gvalues3[k] = gRK3(nu, theta[k])
    gvalues4[k] = gRK4(nu, theta[k])

print("CFL = ", nu)
print("RK1, Max amp factor = ", np.max(np.abs(gvalues1)))
print("RK2, Max amp factor = ", np.max(np.abs(gvalues2)))
print("RK3, Max amp factor = ", np.max(np.abs(gvalues3)))
print("RK4, Max amp factor = ", np.max(np.abs(gvalues4)))

plt.plot(theta, np.abs(gvalues1),label='RK1')
plt.plot(theta, np.abs(gvalues2),label='RK2')
plt.plot(theta, np.abs(gvalues3),label='RK3')
plt.plot(theta, np.abs(gvalues4),label='RK4')
plt.xlabel('Wave number')
plt.ylabel('Amp factor')
plt.legend()
plt.show()
