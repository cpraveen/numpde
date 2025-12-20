'''
Solve u_t + a u_x = 0 in (0,1) with dirichlet BC at x=0
Lax-Wendroff scheme in interior, FTBS at x=1
'''
import numpy as np
import matplotlib.pyplot as plt
import argparse

def f(x):
    return np.sin(2*np.pi*x)

def solve(N, cfl, Tf):
    xmin, xmax = 0.0, 1.0
    a          = 1.0

    h = (xmax - xmin)/N
    dt= cfl * h / np.abs(a)
    nu= a * dt / h

    x = np.linspace(xmin, xmax, N+1)
    u = f(x); unew = np.empty_like(u)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    line1, = ax.plot(x, u, 'r')
    line2, = ax.plot(x, u, 'b')
    ax.set_xlabel('x'); ax.set_ylabel('u')
    plt.grid(True)
    plt.legend(('Numerical','Exact'))
    plt.title('N='+str(N)+', CFL='+str(cfl))
    plt.draw(); plt.pause(0.1)
    wait = input("Press enter to continue ")

    t, it = 0.0, 0
    while t < Tf:
        t += dt; it += 1
        # Left boundary point
        unew[0] = f(xmin-a*t)
        # Interior points, use Lax-Wendroff
        unew[1:-1] = u[1:-1] - 0.5*nu*(u[2:]-u[0:-2]) \
                     + 0.5*nu**2*(u[0:-2]-2*u[1:-1]+u[2:])
        # Last point, use FTBS
        unew[-1] = (1-nu)*u[-1] + nu*u[-2]
        u[:] = unew

        line1.set_ydata(u)
        line2.set_ydata(f(x-a*t))
        plt.draw(); plt.pause(0.1)
    plt.show()

# Get arguments
parser = argparse.ArgumentParser()
parser.add_argument('-N', type=int, help='Number of cells', default=100)
parser.add_argument('-cfl', type=float, help='CFL number', default=0.9)
parser.add_argument('-Tf', type=float, help='Final time', default=1.0)
args = parser.parse_args()

# Run the solver
solve(args.N, args.cfl, args.Tf)
