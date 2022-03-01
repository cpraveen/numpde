"""
Solve diffusion eqn using cell-centered fvm
  -u'' = f in (0,1)
   u(0) = u(1) = 0
with f(x) = sin(pi*x)
"""
from numpy import zeros,linspace,sin,pi,empty_like,log
from numpy.random import rand
import matplotlib.pyplot as plt
import argparse
from prettytable import PrettyTable
from tdma import *

xmin, xmax = 0.0, 1.0   # domain
uexact = lambda x: sin(pi*x)
f   = lambda x: pi**2 * sin(pi*x)

# xf: location of cell faces
def solve(xf,f):
    N  = len(xf) - 1             # number of cells
    xc = 0.5*(xf[0:-1] + xf[1:]) # cell centers
    hc = xf[1:] - xf[0:-1]       # cell lenghts
    hf = empty_like(xf)          # h_{j+1/2}
    hf[0]    = 0.5 * hc[0]
    hf[1:-1] = 0.5*(hc[0:-1] + hc[1:])
    hf[-1]   = 0.5 * hc[-1]

    # Three diagonals
    A,B,C = zeros(N),zeros(N),zeros(N)

    A[0], C[0] = (1/hf[0] + 1/hf[1]), -1/hf[1]
    for i in range(1,N-1):
        B[i],A[i],C[i] = -1/hf[i], (1/hf[i] + 1/hf[i+1]), -1/hf[i+1]
    B[-1], A[-1] = -1/hf[-2], (1/hf[-2] + 1/hf[-1])

    # Right hand side
    F = hc * f(xc)

    # Solve
    u = zeros(N+2)
    u[1:-1] = tdma(A,B,C,F)
    x = zeros(N+2)
    x[0]    = xmin
    x[1:-1] = xc
    x[-1]   = xmax

    return x, u, hc

#------------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('-grid', choices=('uniform','random'), 
                    help='Grid type', required=True)
parser.add_argument('-N', type=int,
                    help='Init no. cells', default=40)
parser.add_argument('-delta', type=float,
                    help='Grid randomness', default=0.1)
args = parser.parse_args()


N     = args.N      # no. of cells
delta = args.delta  # random grid

# Make initial mesh
if args.grid == 'uniform':
    xf = linspace(xmin,xmax,N+1)
else:
    xf = linspace(xmin,xmax,N+1)
    h = xf[1] - xf[0]
    xf[1:-1] += delta * h * (2 * rand(N-1) - 1)

Nlist, Errlist = [], []
table = PrettyTable(['N','hmax','hmax/hmin','err','rate'])
for i in range(5):
    x, u, hc = solve(xf,f)

    # Plot initial mesh solution
    if i == 0:
        plt.figure()
        plt.plot(x,u,'o',label="FVM")
        plt.plot(x,uexact(x),label="Exact")
        plt.legend()
        plt.title('N = '+str(N)+', Grid = '+args.grid)

    # Exact solution to estimate error
    uex = uexact(x)
    hmin = hc.min(); hmax = hc.max(); err = abs(uex-u).max()
    hratio = hmax / hmin
    Nlist.append(N); Errlist.append(err)
    if i == 0:
        table.add_row(['%3d' % N, '%.4e' % hmax, '%.4f' % hratio,
                       '%.4e' % err, '---'])
    else:
        rate = log(Errlist[-2]/Errlist[-1]) / log(Nlist[-1]/Nlist[-2])
        table.add_row(['%3d' % N, '%.4e' % hmax, '%.4f' % hratio,
                       '%.4e' % err, '%.2f' % rate])

    # Refine grid
    N   = 2*N
    tmp = zeros(N + 1)
    tmp[0::2] = xf
    tmp[1::2] = 0.5*(xf[0:-1] + xf[1:])
    xf = tmp

print(table)
plt.figure()
plt.loglog(Nlist,Errlist,'o-')
plt.title('Grid type = '+args.grid)
plt.xlabel('Number of cells')
plt.ylabel('Maximum error')
plt.show()
