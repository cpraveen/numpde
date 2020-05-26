'''
Solve u_t = mu * u_xx in (0,1) with zero dirichlet BC
'''
import numpy as np
import matplotlib.pyplot as plt
import argparse
from ic import *
from tdma import *

xmin, xmax = 0.0, 1.0
mu = 1.0

def solve(scheme, N, lam, Tf, uexact):
    h = (xmax - xmin)/(N-1)
    dt= lam * h**2 / mu

    x = np.linspace(xmin, xmax, N)
    u = uexact(x,0.0)
    u[0], u[-1] = 0.0, 0.0

    fig = plt.figure()
    ax = fig.add_subplot(111)
    line1, = ax.plot(x, u, 'r')
    line2, = ax.plot(x, u, 'b')
    ax.set_xlabel('x'); ax.set_ylabel('u')
    plt.legend(('Exact','Numerical'))
    plt.title('Scheme='+scheme+', N='+str(N)+', $\lambda$='+str(lam))
    plt.draw(); plt.pause(0.1)
    input("Press enter to continue ")

    # Construct tridiagonal matrix and compute its LU decomposition
    if scheme == 'btcs':
        a = (1.0+2.0*lam)*np.ones(N-2)
        b = -lam*np.ones(N-2)
        c = -lam*np.ones(N-2)
        tdma1 (a, b, c)
    elif scheme == 'cn':
        a = (1.0+lam)*np.ones(N-2)
        b = -0.5*lam*np.ones(N-2)
        c = -0.5*lam*np.ones(N-2)
        tdma1 (a, b, c)

    t, it = 0.0, 0
    while t < Tf:
        if scheme == 'ftcs':
            u[1:-1] = lam*u[0:-2] + (1-2*lam)*u[1:-1] + lam*u[2:]
        elif scheme == 'btcs':
            u[1:-1] = tdma2(a, b, c, u[1:-1])
        elif scheme == 'cn':
            rhs = 0.5*lam*u[0:-2] + (1-lam)*u[1:-1] + 0.5*lam*u[2:]
            u[1:-1] = tdma2(a, b, c, rhs)
        t += dt; it += 1
        print("t = ", t)

        line1.set_ydata(uexact(x,t))
        line2.set_ydata(u)
        plt.draw(); plt.pause(0.1)
    plt.show()

# Get arguments
parser = argparse.ArgumentParser()
parser.add_argument('-scheme', choices=('ftcs','btcs','cn'),
                    help='Scheme', required=True)
parser.add_argument('-N', type=int, help='Number of cells', default=51)
parser.add_argument('-lam', type=float, help='Lambda', default=0.49)
parser.add_argument('-Tf', type=float, help='Final time', default=0.025)
parser.add_argument('-ic', choices=('sine1','sine2','const','tri','square'),
                    help='Init cond',   default='sine1')
args = parser.parse_args()

# Select the initial condition function
if args.ic == 'sine1':
    ic = sine1
elif args.ic == 'sine2':
    ic = sine2
elif args.ic == 'const':
    ic = const
elif args.ic == 'square':
    ic = square
elif args.ic == 'tri':
    ic = tri

# Run the solver
solve(args.scheme, args.N, args.lam, args.Tf, ic)
