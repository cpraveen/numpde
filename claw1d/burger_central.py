"""
Solve burger with periodic bc
Central finite difference in space
RK4 for time integration for du/dt = r(u)
"""
import numpy as np
import matplotlib.pyplot as plt
import argparse

def uinit(x):
    return 0.2 * np.sin(x)

def flux(u):
    return 0.5 * u**2

# Returns -f_x using central difference
def rhs_cent(h, u):
    r = np.empty_like(u)
    f = flux(u)
    r[0]    = -(0.5/h)*(f[1]  - f[-2])    # first point
    r[1:-1] = -(0.5/h)*(f[2:] - f[0:-2]) # second to last but one
    r[-1]   = r[0]                       # last is same as first
    return r

# Energy conserving central flux
def flux_ec(ul, ur):
    return (1.0/6.0) * (ul**2 + ul * ur + ur**2)

# Returns -f_x using central difference
def rhs_ec(h, u):
    uim1 = np.roll(u,1)
    ui   = u
    f    = flux_ec(uim1, ui)
    r       = np.empty_like(u)
    r[0:-1] = -(1.0/h)*(f[1:]  - f[0:-1]) # first point
    r[-1]   = r[0]                        # last is same as first
    return r

# Grid has N+1 points with x[0] = x[N]
def solve(N, cfl, Tf, flux):
    if flux == 'cent':
        rhs = rhs_cent
    else:
        rhs = rhs_ec

    xmin, xmax = 0.0, 2.0 * np.pi
    h = (xmax - xmin)/N
    x = np.linspace(xmin, xmax, N+1)
    u = uinit(x)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    line1, = ax.plot(x, u, 'r')
    line2, = ax.plot(x, u, 'b')
    ax.set_xlabel('x'); ax.set_ylabel('u'); plt.grid(True)
    plt.legend(('Solution','Initial'))
    plt.title('N='+str(N)+', CFL='+str(cfl)+', Flux='+flux)
    plt.draw(); plt.pause(0.1)
    wait = input("Press enter to continue ")

    t, it = 0.0, 0
    tdata, edata = [], []
    while t < Tf:
        dt= cfl * h / np.abs(u).max()
        k0 = rhs(h, u)
        k1 = rhs(h, u+0.5*dt*k0)
        k2 = rhs(h, u+0.5*dt*k1)
        k3 = rhs(h, u+dt*k2)
        u += (dt/6)*(k0 + 2*k1 + 2*k2 + k3)
        t += dt; it += 1
        energy = h * np.sum(u[0:-1]**2)
        tdata.append(t); edata.append(energy)
        line1.set_ydata(u)
        plt.ylim(u.min(), u.max())
        plt.draw(); plt.pause(0.1)
    tdata, edata = np.array(tdata), np.array(edata)
    edata = (edata - edata[0])/edata[0]
    plt.figure()
    plt.plot(tdata, edata)
    plt.xlabel('Time, t'); plt.ylabel('Energy, (E(t) - E(0))/E(0)')
    plt.show()

# Get arguments
parser = argparse.ArgumentParser()
parser.add_argument('-N', type=int, help='Number of cells', default=100)
parser.add_argument('-cfl', type=float, help='CFL number', default=0.9)
parser.add_argument('-Tf', type=float, help='Final time', default=10.0)
parser.add_argument('-flux', choices=('cent','ec'), help='Numerical flux', 
                    default='ec')
args = parser.parse_args()

# Run the solver
solve(args.N, args.cfl, args.Tf, args.flux)
