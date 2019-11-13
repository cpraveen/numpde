'''
Solve u_t = u_xx in (0,2*pi) with periodic BC using FTCS scheme
'''
import numpy as np
import matplotlib.pyplot as plt
import argparse

xmin, xmax = 0.0, 2.0*np.pi

def uexact(x,t):
    return 1.0 + np.exp(-t) * np.sin(x) + np.exp(-100*t) * np.sin(10*x)

# N+1 points with h = 2*pi/(N+1), x[0] = 0, x[N] = 2*pi-h
def solve(N, lam, Tf):
    h = (xmax - xmin)/(N+1)
    dt= lam * h**2

    x = xmin + h*np.arange(0,N+1)
    u = uexact(x,0.0);
    unew = np.empty_like(u)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    line1, = ax.plot(x, u, 'r')
    line2, = ax.plot(x, u, 'b')
    ax.set_xlabel('x'); ax.set_ylabel('u')
    plt.legend(('Exact','Numerical'))
    plt.title('N='+str(N)+', $\lambda$='+str(lam))
    plt.draw(); plt.pause(0.1)
    wait = input("Press enter to continue ")

    t, it = 0.0, 0
    while t < Tf:
        if t+dt > Tf:
            dt = Tf-t
            lam = dt/h**2
        # first point
        unew[0] = lam*u[-1] + (1-2*lam)*u[0] + lam*u[1]
        # interior points
        unew[1:-1] = lam*u[0:-2] + (1-2*lam)*u[1:-1] + lam*u[2:]
        # last point
        unew[-1] = lam*u[-2] + (1-2*lam)*u[-1] + lam*u[0]
        t += dt; it += 1
        print("t = ", t)

        line1.set_ydata(uexact(x,t))
        line2.set_ydata(u)
        plt.draw(); plt.pause(0.1)

        u[:] = unew

    plt.figure()
    plt.plot(x,uexact(x,t)-u)
    plt.ylabel('Error')
    plt.xlabel('x')
    plt.title('Error at t='+str(t))
    plt.show()

# Get arguments
parser = argparse.ArgumentParser()
parser.add_argument('-N', type=int, help='Number of cells', default=100)
parser.add_argument('-lam', type=float, help='Lambda', default=0.5)
parser.add_argument('-Tf', type=float, help='Final time', default=0.01)
args = parser.parse_args()

# Run the solver
solve(args.N, args.lam, args.Tf)
