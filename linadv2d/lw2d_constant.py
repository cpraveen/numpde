import numpy as np
import matplotlib.pyplot as plt

nx, ny = 100, 100      # grid size
xmin, xmax = -1.0, 1.0 # domain size
ymin, ymax = -1.0, 1.0
u, v = 1.0, 1.0        # advection velocity
Tf = 2.0               # one time period

# make grid
dx, dy = (xmax - xmin)/nx, (ymax - ymin)/ny
x = np.linspace(xmin+0.5*dx, xmax-0.5*dx, nx)
y = np.linspace(ymin+0.5*dy, ymax-0.5*dy, ny)
Y,X = np.meshgrid(y,x)

# set initial condition
q = np.sin(2*np.pi*X) * np.sin(2*np.pi*Y)

# save copy of ic, we use it to compute error norm
q0 = q.copy()

# plot initial condition
plt.contour(X,Y,q)
plt.xlabel('x'); plt.ylabel('y'); plt.axis('equal')
plt.title('Initial condition')
plt.draw(); plt.pause(1)
wait = input("Press enter to continue ")

# dt from cfl condition
dt = 0.72/(np.abs(u)/dx + np.abs(v)/dy)
print('Grid size nx, ny = ', nx, ny)
print('Grid size dx, dy = ', dx, dy)
print('Time step        = ', dt)

s1, s2 = u*dt/dx, v*dt/dy

# Time loop
t, it = 0.0, 0
while t < Tf:
    if t+dt > Tf: # adjust dt so we reach Tf exactly
        dt = Tf - t
    # generate stencil values using periodicity
    qim1j = np.roll(q,(1,0),(0,1))
    qip1j = np.roll(q,(-1,0),(0,1))
    qijm1 = np.roll(q,(0,1),(0,1))
    qijp1 = np.roll(q,(0,-1),(0,1))
    qim1jm1 = np.roll(q,(1,1),(0,1))
    qim1jp1 = np.roll(q,(1,-1),(0,1))
    qip1jm1 = np.roll(q,(-1,1),(0,1))
    qip1jp1 = np.roll(q,(-1,-1),(0,1))
    # LW scheme
    q = q - 0.5 * s1 * (qip1j - qim1j) - 0.5 * s2 * (qijp1 - qijm1) \
          + 0.5 * s1**2 * (qim1j - 2.0*q + qip1j) \
          + 0.25 * s1 * s2 * (qip1jp1 - qip1jm1 - qim1jp1 + qim1jm1) \
          + 0.5 * s2**2 * (qijm1 - 2.0*q + qijp1)
    t += dt; it += 1
    print('it,t = ',it,t)
    if it%5 == 0 or np.abs(t-Tf) < 1.0e-12:
        plt.clf()
        plt.contour(X,Y,q)
        plt.title('t = '+str(t))
        plt.xlabel('x'); plt.ylabel('y'); plt.axis('equal')
        plt.draw(); plt.pause(0.1)

# Compute error norm: initial condition is exact solution
l1_err = np.sum(np.abs(q-q0)) / (nx*ny)
l2_err = np.sqrt(np.sum((q-q0)**2) / (nx*ny))
li_err = np.abs(q-q0).max()
print('dx,dy,l1,l2,linf error = %10.4e %10.4e %10.4e %10.4e %10.4e' % 
      (dx,dy,l1_err,l2_err,li_err))
plt.show()
