"""
Solve u_t + u_x = 0 with periodic bc
Finite volume scheme with different reconstruction schemes
like first order, minmod, weno5, mp5. Time integration is first order
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

def minmod2(a,b):
    if a*b > 0:
        return np.sign(a) * np.min(np.abs([a,b]))
    else:
        return 0.0

def minmod3(a,b,c):
    if a*b > 0 and b*c > 0:
        return np.sign(a) * np.min(np.abs([a,b,c]))
    else:
        return 0.0

def minmod4(a,b,c,d):
    if a*b > 0 and b*c > 0 and c*d > 0 and d*a > 0:
        return np.sign(a) * np.min(np.abs([a,b,c,d]))
    else:
        return 0.0

def median(a,b,c):
    return a + minmod2(b-a, c-a)

def vanleer(a,b):
    if a*b <= 0.0:
        return 0.0
    else:
        return 2.0*a*b/(a + b)

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

# Scheme of Suresh and Huynh
def mp5(um2,um1,u0,up1,up2):
    eps, alpha = 1.0e-13, 4.0
    u = (2.0*um2 - 13.0*um1 + 47.0*u0 + 27.0*up1 - 3.0*up2)/60.0
    ump = u0 + minmod2(up1-u0, alpha*(u0-um1))
    if (u - u0)*(u - ump) < eps:
        return u
    d0 = um1 + up1 - 2.0*u0
    dm1= um2 + u0  - 2.0*um1
    dp1= u0  + up2 - 2.0*up1
    dlm4 = minmod4(4*dm1 - d0, 4*d0-dm1, dm1, d0)
    drm4 = minmod4(4*d0 - dp1, 4*dp1-d0, d0,  dp1)
    uul = u0 + alpha*(u0 - um1)
    uav = 0.5*(u0 + up1)
    umd = uav - 0.5*drm4
    ulc = u0 + 0.5*(u0 - um1) + (4.0/3.0)*dlm4
    umin = np.max([np.min([u0,up1,umd]), np.min([u0,uul,ulc])])
    umax = np.min([np.max([u0,up1,umd]), np.max([u0,uul,ulc])])
    u = median(u, umin, umax)
    return u


# Compute left state at interface b/w uj and ujp1
def reconstruct(ujm2, ujm1, uj, ujp1, ujp2):
    if scheme==0: # first order
        return uj
    elif scheme==1: # minmod
        return uj + 0.5 * minmod3(beta*(uj-ujm1), 0.5*(ujp1-ujm1), beta*(ujp1-uj))
    elif scheme==2: # weno5
        return weno5(ujm2, ujm1, uj, ujp1, ujp2)
    elif scheme==3: # mp5
        return mp5(ujm2, ujm1, uj, ujp1, ujp2)
    elif scheme==4: # van leer
        return uj + 0.5 * vanleer(uj-ujm1, ujp1-uj)

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
def solve(N, cfl, rscheme, Tf, xmin, xmax, uinit, nrk):
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
    plt.axis([xmin, xmax, u.min()-0.1, u.max()+0.1])
    plt.grid(True); plt.draw(); plt.pause(0.1)
    wait = input("Press enter to continue ")

    t, it = 0.0, 0
    while t < Tf:
        uold[:] = u
        for rk in range(nrk):
            res = compute_residual(u)
            u = ark[rk]*uold + (1.0-ark[rk])*(u - (dt/h)*res)
        t += dt; it += 1
        line1.set_ydata(u)
        line2.set_ydata(uinit(x-a*t))
        plt.draw(); plt.pause(0.1)
    plt.show()

# Get arguments
parser = argparse.ArgumentParser()
parser.add_argument('-N', type=int, help='Number of cells', default=100)
parser.add_argument('-cfl', type=float, help='CFL number', default=0.9)
parser.add_argument('-scheme', choices=('FO','MMOD','VL','WENO5','MP5'),
                    help='Scheme', default='FO')
parser.add_argument('-Tf', type=float, help='Final time', default=1.0)
parser.add_argument('-ic', choices=('smooth','hat','mult'), help='Init cond', default='smooth')
args = parser.parse_args()

if args.scheme=="FO":
    scheme, nrk = 0, 1
elif args.scheme=="MMOD":
    scheme, nrk = 1, 3
elif args.scheme=="WENO5":
    scheme, nrk = 2, 3
elif args.scheme=="MP5":
    scheme, nrk = 3, 3
elif args.scheme == "VL":
    scheme, nrk = 4, 3

# Run the solver
if args.ic == "smooth":
    xmin, xmax = 0.0, 1.0
    solve(args.N, args.cfl, args.scheme, args.Tf, xmin, xmax, smooth, nrk)
elif args.ic == "hat":
    xmin, xmax = 0.0, 1.0
    solve(args.N, args.cfl, args.scheme, args.Tf, xmin, xmax, hat, nrk)
elif args.ic == "mult":
    xmin, xmax = -1.0, 1.0
    solve(args.N, args.cfl, args.scheme, args.Tf, xmin, xmax, mult, nrk)
