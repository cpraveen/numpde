'''
Solve u_t = mu * u_xx in (0,1) with dirichlet BC
using FTCS scheme
'''
import numpy as np
import matplotlib.pyplot as plt
import argparse

xmin, xmax = 0.0, 1.0
mu = 1.0

def uexact(x,t):
    return np.exp(-mu*np.pi**2*t) * np.sin(np.pi*x)

def solve(N, lam, Tf):
    h = (xmax - xmin)/(N-1)
    dt= lam * h**2 / mu

    x = np.linspace(xmin, xmax, N)
    u = uexact(x,0.0);

    fig = plt.figure()
    ax = fig.add_subplot(111)
    line1, = ax.plot(x, u, 'r')
    line2, = ax.plot(x, u, 'b')
    ax.set_xlabel('x'); ax.set_ylabel('u')
    plt.legend(('Exact','Numerical'))
    plt.title('N='+str(N)+', $\lambda$='+str(lam))
    plt.draw(); plt.pause(0.1)
    wait = raw_input("Press enter to continue ")

    t, it = 0.0, 0
    while t < Tf:
        u[1:-1] = lam*u[0:-2] + (1-2*lam)*u[1:-1] + lam*u[2:]
        t += dt; it += 1
        print "t = ", t

        line1.set_ydata(uexact(x,t))
        line2.set_ydata(u)
        plt.draw(); plt.pause(0.1)
    plt.show()

# Get arguments
parser = argparse.ArgumentParser()
parser.add_argument('-N', type=int, help='Number of cells', default=51)
parser.add_argument('-lam', type=float, help='Lambda', default=0.49)
parser.add_argument('-Tf', type=float, help='Final time', default=0.025)
args = parser.parse_args()

# Run the solver
solve(args.N, args.lam, args.Tf)
