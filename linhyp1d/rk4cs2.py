"""
Solve u_t + u_x = 0 with periodic bc
Central finite difference in space
RK4 for time integration for du/dt = r(u)
"""
import numpy as np
import matplotlib.pyplot as plt
import argparse
from ic import *

# Returns -u_x using central difference
def rhs(h, u):
    r = np.empty_like(u)
    r[0]    = -(0.5/h)*(u[1] - u[-2])    # first point
    r[1:-1] = -(0.5/h)*(u[2:] - u[0:-2]) # second to last but one
    r[-1]   = r[0]                       # last is same as first
    return r

def solve(N, cfl, Tf, uinit):
    xmin, xmax = 0.0, 1.0
    a          = 1.0

    h = (xmax - xmin)/N
    dt= cfl * h / np.abs(a)

    x = np.linspace(xmin, xmax, N+1)
    u = uinit(x)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    line1, = ax.plot(x, u, 'r')
    line2, = ax.plot(x, u, 'b')
    ax.set_xlabel('x'); ax.set_ylabel('u')
    plt.legend(('Numerical','Exact'))
    plt.title('N='+str(N)+', CFL='+str(cfl)+', Scheme=RK4+CS')
    plt.draw(); plt.pause(0.1)
    wait = raw_input("Press enter to continue ")

    t, it = 0.0, 0
    while t < Tf:
        k0 = a * rhs(h, u)
        k1 = a * rhs(h, u+0.5*dt*k0)
        k2 = a * rhs(h, u+0.5*dt*k1)
        k3 = a * rhs(h, u+dt*k2)
        u += (dt/6)*(k0 + 2*k1 + 2*k2 + k3)
        t += dt; it += 1
        line1.set_ydata(u)
        line2.set_ydata(uinit(x-a*t))
        plt.draw(); plt.pause(0.1)
    plt.show()

# Get arguments
parser = argparse.ArgumentParser()
parser.add_argument('-N', type=int, help='Number of cells', default=100)
parser.add_argument('-cfl', type=float, help='CFL number', default=0.9)
parser.add_argument('-Tf', type=float, help='Final time', default=1.0)
parser.add_argument('-ic', choices=('smooth','hat'), help='Init cond', default='smooth')
args = parser.parse_args()

# Run the solver
if args.ic == "smooth":
    solve(args.N, args.cfl, args.Tf, smooth)
else:
    solve(args.N, args.cfl, args.Tf, hat)
