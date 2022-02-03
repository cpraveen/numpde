'''
Solve variable advection with fourier spectral method and leap-frog in time.
Same problem as in p06 in Trefethen: Spectral methods in Matlab
See : https://nbviewer.org/github/cpraveen/chebpy/blob/master/p06.ipynb
'''
from numpy import pi,linspace,sin,exp,round,zeros,arange,real,abs
from numpy.fft import fft,ifft
import matplotlib.pyplot as plt

# Initial condition
f = lambda x: exp(-100.0*(x-1.0)**2.0)

# Some parameters
cfl = 0.1  # cfl number
N   = 128  # number of grid points

h = 2*pi/N
x = h*arange(1,N+1);
c = 0.2 + sin(x-1.0)**2.0
dt = cfl * h / abs(c).max()
v = f(x); vold = f(x-c*dt)

fig = plt.figure()
ax = fig.add_subplot(111)
line, = ax.plot(x, v, 'r')
ax.set_xlabel('x'); ax.set_ylabel('v')
plt.grid(True)
plt.title('N='+str(N))
plt.draw(); plt.pause(0.1)
wait = input("Press enter to continue ")

# Time-stepping by leap-frog formula
t, tmax = 0.0, 8.0
niter= int(tmax / dt)
tdata = []; tdata.append(0.0)
for i in range(niter):
    t = t + dt
    v_hat = fft(v)
    w_hat = 1j*zeros(N)
    w_hat[0:N//2] = 1j*arange(0,N//2)
    w_hat[N//2+1:] = 1j*arange(-N//2+1,0,1)
    w_hat = w_hat * v_hat
    w = real(ifft(w_hat))
    vnew = vold - 2.0*dt*c*w
    vold = v; v = vnew;
    if i % 10 == 0:
        line.set_ydata(v)
        plt.draw(); plt.pause(0.1)

plt.show()
