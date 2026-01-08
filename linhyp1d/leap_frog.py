"""
Solve u_t + u_x = 0 with periodic bc
"""
import numpy as np
import matplotlib.pyplot as plt
import argparse
from ic import *
from errornorm import *

def solve(N, cfl, Tf, uinit):
    xmin, xmax = 0.0, 1.0
    a          = 1.0

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
    plt.title('N='+str(N)+', CFL='+str(cfl)+', Scheme=Leap-Frog')
    plt.draw(); plt.pause(0.1)
    wait = input("Press enter to continue ")

    t, it, times, error, = 0.0, 0, [], []

    # First step: FTCS
    u1 = np.empty_like(u)
    u1[0]    = u[0] + 0.5*nu*(u[-2] - u[1])
    u1[1:-1] = u[1:-1] + 0.5*nu*(u[0:-2] - u[2:])
    u1[-1]   = u1[0]

    t += dt; it += 1
    times.append(t)
    error.append(errornorm(h, u1, uinit(x-a*t)))

    # Now use Leap-Frog: u -> u^(n-1), u1 -> u^(n), u2 -> u^(n+1)
    u2 = np.empty_like(u)
    while t < Tf:
        u2[0]    = u[0] + nu*(u1[-2] - u1[1])
        u2[1:-1] = u[1:-1] + nu*(u1[0:-2] - u1[2:])
        u2[-1]   = u2[0]
        t += dt; it += 1
        line1.set_ydata(u2)
        line2.set_ydata(uinit(x-a*t))
        plt.draw(); plt.pause(0.1)
        u[:]  = u1
        u1[:] = u2
        times.append(t)
        error.append(errornorm(h, u2, uinit(x-a*t)))

    np.savetxt('error_LEAP.txt',np.column_stack([times,error]))
    plt.figure()
    plt.plot(times, error)
    plt.xlabel('t'); plt.ylabel('Error norm')
    plt.grid(True)
    plt.show()

# Get arguments
parser = argparse.ArgumentParser()
parser.add_argument('-N', type=int, help='Number of cells', default=100)
parser.add_argument('-cfl', type=float, help='CFL number', default=0.98)
parser.add_argument('-Tf', type=float, help='Final time', default=1.0)
parser.add_argument('-ic', choices=('smooth','hat','sine'), help='Init cond', default='smooth')
args = parser.parse_args()

# Run the solver
if args.ic == "smooth":
    solve(args.N, args.cfl, args.Tf, smooth)
elif args.ic == "hat":
    solve(args.N, args.cfl, args.Tf, hat)
elif args.ic == "sine":
    solve(args.N, args.cfl, args.Tf, sine)
