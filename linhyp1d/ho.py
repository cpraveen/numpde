"""
Solve u_t + u_x = 0 with periodic bc
Finite volume scheme with different reconstruction schemes
like first order, minmod, weno5. Time integration is first order
euler or ssprk3
"""
import numpy as np
import matplotlib.pyplot as plt
import argparse
from ic import *

# Coefficients of 3-stage SSPRK
ark = [0.0, 3.0/4.0, 1.0/3.0]
a = 1.0
beta = 2.0 # used in minmod scheme

def minmod(a,b,c):
    if a*b > 0 and b*c > 0:
        return np.sign(a) * np.min(np.abs([a,b,c]))
    else:
        return 0.0

# 5th order weno reconstruction. Gives left state at
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

# Compute left state at interface b/w uj and ujp1
def reconstruct(ujm2, ujm1, uj, ujp1, ujp2):
    if scheme==0: # first order
        return uj
    elif scheme==1: # minmod
        return uj + 0.5 * minmod(beta*(uj-ujm1), 0.5*(ujp1-ujm1), beta*(ujp1-uj))
    elif scheme==2: # weno5
        return weno5(ujm2, ujm1, uj, ujp1, ujp2)

# Compute finite volume residual R in eqn h*du/dt + R = 0
def compute_residual(u):
    n = len(u)
    res = np.zeros(n)

    # There are n+1 faces. Compute flux across 0,...,n-1 face
    # Flux across n'th face = flux across 0'th face, so no need
    # to compute it again. Note that res[-1] is residual of last cell
    # etc.
    for f in range(n-1):
        ul = reconstruct(u[f-3], u[f-2], u[f-1], u[f], u[f+1])
        flux = a * ul
        res[f-1] += flux
        res[f]   -= flux

    # penultimate face
    f = n-1
    ul = reconstruct(u[f-3], u[f-2], u[f-1], u[f], u[0])
    flux = a * ul
    res[f-1] += flux
    res[f]   -= flux

    return res

# Solve the problem
def solve(N, cfl, rscheme, Tf, uinit, nrk):
    xmin, xmax = 0.0, 1.0

    h = (xmax - xmin)/N
    dt= cfl * h / np.abs(a)

    x = np.linspace(xmin+0.5*h, xmax-0.5*h, N)
    u = uinit(x)
    uold = np.empty_like(u)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    line1, = ax.plot(x, u, 'ro')
    line2, = ax.plot(x, u, 'b')
    ax.set_xlabel('x'); ax.set_ylabel('u')
    plt.legend(('Numerical','Exact'))
    plt.title('N='+str(N)+', CFL='+str(cfl)+', Scheme='+rscheme)
    plt.axis([0.0, 1.0, -0.1, 1.1])
    plt.draw(); plt.pause(0.1)
    wait = raw_input("Press enter to continue ")

    t, it = 0.0, 0
    while t < Tf:
        uold[:] = u
        for rk in range(nrk):
            res = compute_residual(u)
            u = ark[rk]*uold + (1-ark[rk])*(u - (dt/h)*res)
        t += dt; it += 1
        line1.set_ydata(u)
        line2.set_ydata(uinit(x-a*t))
        plt.draw(); plt.pause(0.1)
    plt.show()

# Get arguments
parser = argparse.ArgumentParser()
parser.add_argument('-N', type=int, help='Number of cells', default=100)
parser.add_argument('-cfl', type=float, help='CFL number', default=0.9)
parser.add_argument('-scheme', choices=('FO','MMOD','WENO'), help='Scheme', default='FO')
parser.add_argument('-Tf', type=float, help='Final time', default=1.0)
parser.add_argument('-ic', choices=('smooth','hat'), help='Init cond', default='smooth')
args = parser.parse_args()

if args.scheme=="FO":
    scheme, nrk = 0, 1
elif args.scheme=="MMOD":
    scheme, nrk = 1, 3
elif args.scheme=="WENO":
    scheme, nrk = 2, 3

# Run the solver
if args.ic == "smooth":
    solve(args.N, args.cfl, args.scheme, args.Tf, smooth, nrk)
else:
    solve(args.N, args.cfl, args.scheme, args.Tf, hat, nrk)
