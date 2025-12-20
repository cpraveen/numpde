'''
Plots amplification factor for forward Euler + second order backward difference
scheme for different CFL numbers. All the curves have portions which lie above
unity and the scheme is not stable for any CFL.
'''
import numpy as np
import matplotlib.pyplot as plt

s = [0.01, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8]
kh = np.linspace(-np.pi, np.pi, 500)

plt.figure()
leg = ()
for cfl in s:
    g = 1.0 - 1.5*cfl + 2.0*cfl*np.exp(-1j*kh) - 0.5*cfl*np.exp(-2j*kh)
    g = np.abs(g)
    plt.plot(kh,g)
    leg += ('CFL='+str(cfl),)
    print('cfl, max amp = %6.3f %16.12f' % (cfl, g.max()))

plt.xlabel('Wave number $k h$')
plt.ylabel('Amplification factor $|\gamma_k|$')
plt.legend(leg)
plt.title('Second order upwind with forward Euler')
plt.grid(True)
plt.show()
