import numpy as np
import matplotlib.pyplot as plt

N = 500
zr = np.linspace(-3.0,1.0,N)
zi = np.linspace(-3.0,3.0,N)
zr,zi = np.meshgrid(zr,zi)
z = zr + 1j * zi
G = 1.0 + z + z**2/2 + z**3/6
plt.contour(zr,zi,np.abs(G),levels=[0.9999,1,1.0001],linewidths=2)
plt.plot([0,0],[-1.73,1.73],'sg')

k = np.linspace(0.0, 2*np.pi, 500)
g = -(1.5 - 2 * np.exp(-1j*k) + 0.5 * np.exp(-2j*k))

cfl = 0.75
plt.plot(np.real(cfl*g), np.imag(cfl*g), label='cfl='+str(cfl))

cfl = 0.625
plt.plot(np.real(cfl*g), np.imag(cfl*g), label='cfl='+str(cfl))

cfl = 0.5
plt.plot(np.real(cfl*g), np.imag(cfl*g), label='cfl='+str(cfl))

plt.legend()
plt.grid(True)
plt.title('Second order upwind + RK3')
plt.xlim(-5,1)
plt.ylim(-3,3)
plt.axis('equal')

plt.show()
