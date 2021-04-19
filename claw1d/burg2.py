"""
Solve u_t + f(u)_x = 0  for f(u) = u^2/2
Finite volume scheme
"""
import numpy as np
import matplotlib.pyplot as plt
import argparse
from numfluxes import *

ul, ur = 1.2, 0.4
s      = 0.5*(ul + ur)

def uexact(t,x):
    u = np.zeros(len(x))
    for i in range(len(x)):
        if x[i] < 0.25 + s*t:
            u[i] = ul
        else:
            u[i] = ur
    return u

def solve(N, cfl, scheme, Tf):
    xmin, xmax = 0.0, 1.0

    x = np.linspace(xmin, xmax, N)
    h = (xmax - xmin)/(N-1)
    u = uexact(0.0, x)
    dt= cfl * h / np.max(u)
    lam = dt/h

    fig = plt.figure()
    ax = fig.add_subplot(111)
    line1, = ax.plot(x, u, 'o')
    line2, = ax.plot(x, u, 'r')
    ax.set_xlabel('x'); ax.set_ylabel('u')
    plt.legend(('Numerical','Exact'))
    plt.title('N='+str(N)+', CFL='+str(cfl)+', Scheme='+scheme)
    plt.axis([0.0, 1.0, 0.0, 1.4]); plt.grid(True)
    plt.draw(); plt.pause(0.1)
    wait = input("Press enter to continue ")

    t, it = 0.0, 0
    while t < Tf:
        if scheme=='LF':
            f = flux_lf(lam, u)
        elif scheme=='LLF':
            f = flux_llf(u)
        elif scheme=='LW':
            f = flux_lw(lam, u)
        elif scheme=='ROE':
            f = flux_roe(u)
        elif scheme=='GOD':
            f = flux_god(u)
        else:
            print("Unknown scheme: ", scheme)
            return
        u[1:-1] -= lam * (f[2:-1] - f[1:-2])
        t += dt; it += 1
        line1.set_ydata(u)
        line2.set_ydata(uexact(t,x))
        plt.draw(); plt.pause(0.1)
    fname = 'burg2_'+scheme+'.txt'
    np.savetxt(fname, np.column_stack([x, u, uexact(t,x)]))
    print('Saved file ', fname)
    plt.show()

# Get arguments
parser = argparse.ArgumentParser()
parser.add_argument('-N', type=int, help='Number of cells', default=100)
parser.add_argument('-cfl', type=float, help='CFL number', default=0.9)
parser.add_argument('-scheme', choices=('LF','LLF','LW','ROE','GOD'), 
                    help='Scheme', default='LF')
parser.add_argument('-Tf', type=float, help='Final time', default=0.5)
args = parser.parse_args()

# Run the solver
solve(args.N, args.cfl, args.scheme, args.Tf)
