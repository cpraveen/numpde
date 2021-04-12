import numpy as np
import matplotlib.pyplot as plt

g = lambda cfl, k: (1.0 - cfl*(1.5 - 2 * np.exp(-1j*k) + 0.5 * np.exp(-2j*k)) 
                    + 0.5 * cfl**2 * (np.exp(-2j*k)  - 2 * np.exp(-1j*k) + 1.0))

N = 500
k = np.linspace(0.0, 2*np.pi, N)

plt.plot(np.cos(k), np.sin(k), '--', lw=2, label='Unit circle')

for cfl in [0.5, 1.0, 1.5, 2.0]:
    gg = g(cfl,k)
    plt.plot(np.real(gg), np.imag(gg), label='cfl='+str(cfl))

plt.legend()
plt.grid(True)
plt.title('Beam-Warming method')
#plt.xlim(-5,1)
#plt.ylim(-3,3)
plt.axis('equal')

plt.show()
