"""
Solve u_t + u_x = 0 with periodic bc
"""
import numpy as np
import matplotlib.pyplot as plt
import argparse
from ic import *

# FTBS
def update_ftbs(nu, u):
    unew = np.empty_like(u)
    unew[1:] = (1-nu)*u[1:] + nu*u[0:-1]
    unew[0] = unew[-1]
    return unew

# FTFS
def update_ftfs(nu, u):
    unew = np.empty_like(u)
    unew[0:-1] = (1+nu)*u[0:-1] - nu*u[1:]
    unew[-1] = unew[0]
    return unew

# FTCS
def update_ftcs(nu, u):
    unew = np.empty_like(u)
    unew[0] = u[0] + 0.5*nu*(u[-2] - u[1])
    unew[1:-1] = u[1:-1] + 0.5*nu*(u[0:-2] - u[2:])
    unew[-1] = unew[0]
    return unew

# Lax-Friedrich
def update_lf(nu, u):
    unew = np.empty_like(u)
    unew[0] = 0.5*(u[-2] + u[1]) + 0.5*nu*(u[-2] - u[1])
    unew[1:-1] = 0.5*(u[0:-2] + u[2:]) + 0.5*nu*(u[0:-2] - u[2:])
    unew[-1] = unew[0]
    return unew

# Lax-Wendroff
def update_lw(nu, u):
    unew = np.empty_like(u)
    unew[0] = u[0] - 0.5*nu*(u[1]-u[-2]) + 0.5*nu**2*(u[-2]-2*u[0]+u[1])
    unew[1:-1] = u[1:-1] - 0.5*nu*(u[2:] - u[0:-2]) \
                 + 0.5*nu**2*(u[0:-2] - 2*u[1:-1] + u[2:])
    unew[-1] = unew[0]
    return unew

# Beam-Warming scheme
def update_bw(nu, u):
    unew = np.empty_like(u)
    unew[0] = u[0] - nu*(1.5*u[0] - 2.0*u[-2] + 0.5*u[-3]) \
                   + 0.5*nu**2*(u[-3] - 2.0*u[-2] + u[0])
    unew[1] = u[1] - nu*(1.5*u[1] - 2.0*u[0] + 0.5*u[-2]) \
                   + 0.5*nu**2*(u[-2] - 2.0*u[0] + u[1])
    unew[2:] = u[2:] - nu*(1.5*u[2:] - 2.0*u[1:-1] + 0.5*u[0:-2]) \
                     + 0.5*nu**2*(u[0:-2] - 2.0*u[1:-1] + u[2:])
    return unew

# Solve the problem
def solve(a, N, cfl, scheme, Tf, uinit):
    xmin, xmax = 0.0, 1.0

    h = (xmax - xmin)/N
    dt= cfl * h / np.abs(a)
    nu= a * dt / h

    x = np.linspace(xmin, xmax, N+1)
    u = uinit(x)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    line1, = ax.plot(x, u, 'r')
    line2, = ax.plot(x, u, 'b')
    ax.set_xlabel('x'); ax.set_ylabel('u')
    plt.grid(True)
    plt.legend(('Numerical','Exact'))
    plt.title('N='+str(N)+', CFL='+str(cfl)+', Scheme='+scheme)
    plt.draw(); plt.pause(0.1)
    wait = input("Press enter to continue ")

    t, it = 0.0, 0
    while t < Tf:
        if scheme=='FTBS':
            u = update_ftbs(nu, u)
        elif scheme=='FTFS':
            u = update_ftfs(nu, u)
        elif scheme=='FTCS':
            u = update_ftcs(nu, u)
        elif scheme=='LF':
            u = update_lf(nu, u)
        elif scheme=='LW':
            u = update_lw(nu, u)
        elif scheme=='BW':
            u = update_bw(nu, u)
        else:
            print("Unknown scheme: ", scheme)
            return
        t += dt; it += 1
        line1.set_ydata(u)
        line2.set_ydata(uinit(x-a*t))
        plt.draw(); plt.pause(0.1)
    # Plot final solution with symbols
    plt.figure()
    plt.plot(x,u,'ro',label='Numerical')
    plt.plot(x,uinit(x-a*t),'b-',label='Exact')
    plt.title('N='+str(N)+', CFL='+str(cfl)+', Scheme='+scheme)
    plt.xlabel('x'); plt.ylabel('u')
    plt.legend(); plt.grid(True)
    np.savetxt('sol.txt',np.column_stack([x,u,uinit(x-a*t)]))
    plt.show()

#------------------------------------------------------------------------------
# Get arguments
parser = argparse.ArgumentParser()
parser.add_argument('-N', type=int, help='Number of cells', default=100)
parser.add_argument('-cfl', type=float, help='CFL number', default=0.98)
parser.add_argument('-scheme', choices=('FTBS','FTFS','FTCS','LF','LW','BW'),
                    help='Scheme', default='FTBS')
parser.add_argument('-a', type=float, help='Advection speed', default=1.0)
parser.add_argument('-Tf', type=float, help='Final time', default=1.0)
parser.add_argument('-ic', choices=('smooth','hat','sine'), help='Init cond', default='smooth')
args = parser.parse_args()

# Run the solver
if args.ic == "smooth":
    solve(args.a, args.N, args.cfl, args.scheme, args.Tf, smooth)
elif args.ic == "hat":
    solve(args.a, args.N, args.cfl, args.scheme, args.Tf, hat)
elif args.ic == "sine":
    solve(args.a, args.N, args.cfl, args.scheme, args.Tf, sine)
