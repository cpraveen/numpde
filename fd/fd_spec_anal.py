'''
See https://en.wikipedia.org/wiki/Finite_difference_coefficient
'''
import numpy as np
import matplotlib.pyplot as plt

c2 = [-1.0/2.0, 0.0, 1.0/2.0]
c4 = [1.0/12.0, -2.0/3.0, 0.0, 2.0/3.0, -1.0/12.0]
c6 = [-1.0/60.0, 3.0/20.0, -3.0/4.0, 0.0, 3.0/4.0, -3.0/20.0, 1.0/60.0]

n = 100
w = np.linspace(1.0e-6, np.pi, n)
fd2 = (np.exp(-1j*w)*c2[0] + np.exp(0j*w)*c2[1] + np.exp(1j*w)*c2[2])/w
fd4 = (np.exp(-2j*w)*c4[0] + np.exp(-1j*w)*c4[1] 
     + np.exp(0j*w)*c4[2] + np.exp(1j*w)*c4[3] + np.exp(2j*w)*c4[4])/w
fd6 = (np.exp(-3j*w)*c6[0] + np.exp(-2j*w)*c6[1] + np.exp(-1j*w)*c6[2] 
     + np.exp(0j*w)*c6[3] + np.exp(1j*w)*c6[4] + np.exp(2j*w)*c6[5] 
     + np.exp(3j*w)*c6[6])/w

plt.plot(w, np.abs(fd2), label='2nd order')
plt.plot(w, np.abs(fd4), label='4th order')
plt.plot(w, np.abs(fd6), label='6th order')
plt.plot(w, np.ones(n), '--', label='Ideal')
plt.xlabel('w')
plt.ylabel('$\\varepsilon_d$')
plt.legend()
plt.title('Fourier amplitude of FD schemes')

plt.show()
