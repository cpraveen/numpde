import numpy as np
import matplotlib.pyplot as plt

nx, ny = 100, 100
xmin, xmax = -1.0, 1.0
ymin, ymax = -1.0, 1.0
u, v = 1.0, 1.0
Tf = 2.0

dx, dy = (xmax - xmin)/nx, (ymax - ymin)/ny
x = np.linspace(xmin+0.5*dx, xmax-0.5*dx, nx)
y = np.linspace(ymin+0.5*dy, ymax-0.5*dy, ny)
X,Y = np.meshgrid(x,y)
q = np.sin(2*np.pi*X) * np.sin(2*np.pi*Y)

plt.contour(X,Y,q)
plt.xlabel('x'); plt.ylabel('y'); plt.axis('equal')
plt.title('Initial condition')
plt.draw(); plt.pause(1)
wait = input("Press enter to continue ")

dt = 0.72/(np.abs(u)/dx + np.abs(v)/dy)
print('Grid size nx, ny = ', nx, ny)
print('Grid size dx, dy = ', dx, dy)
print('Time step        = ', dt)

s1, s2 = u*dt/dx, v*dt/dy

t, it = 0.0, 0
while t < Tf:
    if t+dt > Tf:
        dt = Tf - t
    qim1j = np.roll(q,(1,0),(0,1))
    qip1j = np.roll(q,(-1,0),(0,1))
    qijm1 = np.roll(q,(0,1),(0,1))
    qijp1 = np.roll(q,(0,-1),(0,1))
    qim1jm1 = np.roll(q,(1,1),(0,1))
    qim1jp1 = np.roll(q,(1,-1),(0,1))
    qip1jm1 = np.roll(q,(-1,1),(0,1))
    qip1jp1 = np.roll(q,(-1,-1),(0,1))
    q = q - 0.5*s1*(qip1j - qim1j) - 0.5*s2*(qijp1 - qijm1) \
          + 0.5*s1**2*(qim1j - 2.0*q + qip1j) \
          + 0.25*s1*s2*(qip1jp1 - qip1jm1 - qim1jp1 + qim1jm1) \
          + 0.5*s2**2*(qijm1 - 2.0*q + qijp1)
    t += dt; it += 1
    print('it,t = ',it,t)
    if it%5 == 0 or np.abs(t-Tf) < 1.0e-12:
        plt.close()
        plt.contour(X,Y,q)
        plt.title('t = '+str(t))
        plt.xlabel('x'); plt.ylabel('y'); plt.axis('equal')
        plt.draw(); plt.pause(0.1)

plt.show()
