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
    g2 = 1 - 6*cfl + 26*cfl**2 - 8*cfl*(-1+4*cfl)*np.cos(kh) \
            + 2*cfl*(-1+3*cfl)*np.cos(2*kh)
    g2 = np.sqrt(g2)
    plt.plot(kh,g2)
    leg += (str(cfl),)

plt.xlabel('$k h$')
plt.ylabel('$|\gamma_k|$')
plt.legend(leg)
plt.title('Second order upwind with forward Euler')
plt.grid(True)
plt.show()
