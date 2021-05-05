"""
Solve u_t + f(u)_x = 0  for f(u) = u^2/2
Finite volume scheme
"""
import numpy as np
import matplotlib.pyplot as plt
import argparse
from numfluxes import *

def shock(x):
    u = np.zeros(len(x))
    for i in range(len(x)):
        if x[i] < 0.25:
            u[i] = 2.0
        else:
            u[i] = 1.0
    return u

def smooth(x):
    return np.sin(2*np.pi*x)

# Rarefaction without sonic point
def rare1(x):
    u = np.zeros(len(x))
    for i in range(len(x)):
        if x[i] < 0.25:
            u[i] = 0.5
        else:
            u[i] = 2.0
    return u

# Rarefaction with sonic point
def rare(x):
    u = np.zeros(len(x))
    for i in range(len(x)):
        if x[i] < 0.5:
            u[i] = -0.5
        else:
            u[i] = 1.0
    return u

def expo(x):
    return 1.0 + np.exp(-100*(x-0.25)**2)

def slope(x):
    x1, x2 = 0.2, 0.4
    u1, u2 = 1.0, 0.5
    u = np.zeros(len(x))
    for i in range(len(x)):
        if x[i] < x1:
            u[i] = u1
        elif x[i] > x2:
            u[i] = u2
        else:
            u[i] = u1*(x[i] - x2)/(x1 - x2) + u2*(x[i] - x1)/(x2 - x1)
    return u

def solve(N, cfl, scheme, Tf, uinit):
    xmin, xmax = 0.0, 1.0

    x = np.linspace(xmin, xmax, N)
    h = (xmax - xmin)/(N-1)
    u = uinit(x)
    dt= cfl * h / np.max(u)
    lam = dt/h

    fig = plt.figure()
    ax = fig.add_subplot(111)
    line1, = ax.plot(x, u, 'o')
    ax.set_xlabel('x'); ax.set_ylabel('u')
    plt.title('N='+str(N)+', CFL='+str(cfl)+', Scheme='+scheme)
    plt.grid(True); plt.draw(); plt.pause(0.1)
    wait = input("Press enter to continue ")

    t, it = 0.0, 0
    while t < Tf:
        if scheme=='C':
            f = flux_central(lam, u)
        elif scheme=='LF':
            f = flux_lf(lam, u)
        elif scheme=='GLF':
            f = flux_glf(u)
        elif scheme=='LLF':
            f = flux_llf(u)
        elif scheme=='LW':
            f = flux_lw(lam, u)
        elif scheme=='ROE':
            f = flux_roe(u)
        elif scheme=='EROE':
            f = flux_eroe(u)
        elif scheme=='GOD':
            f = flux_god(u)
        else:
            print("Unknown scheme: ", scheme)
            return
        u[1:-1] -= lam * (f[2:-1] - f[1:-2])
        t += dt; it += 1
        line1.set_ydata(u)
        plt.draw(); plt.pause(0.1)
    plt.show()

# Get arguments
parser = argparse.ArgumentParser()
parser.add_argument('-N', type=int, help='Number of cells', default=100)
parser.add_argument('-cfl', type=float, help='CFL number', default=0.9)
parser.add_argument('-scheme',
                    choices=('C','LF','GLF','LLF','LW','ROE','EROE','GOD'),
                    help='Scheme', default='LF')
parser.add_argument('-ic',
                    choices=('smooth','shock','rare1','rare','expo','slope'),
                    help='Initial condition', default='smooth')
parser.add_argument('-Tf', type=float, help='Final time', default=1.0)
args = parser.parse_args()

# Run the solver
if args.ic == "smooth":
    solve(args.N, args.cfl, args.scheme, args.Tf, smooth)
elif args.ic == "expo":
    solve(args.N, args.cfl, args.scheme, args.Tf, expo)
elif args.ic == "shock":
    solve(args.N, args.cfl, args.scheme, args.Tf, shock)
elif args.ic == "rare1":
    solve(args.N, args.cfl, args.scheme, args.Tf, rare1)
elif args.ic == "rare":
    solve(args.N, args.cfl, args.scheme, args.Tf, rare)
elif args.ic == "slope":
    solve(args.N, args.cfl, args.scheme, args.Tf, slope)
