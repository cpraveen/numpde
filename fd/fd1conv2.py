"""
Convergence of finite difference approximation of first derivative when
function is in C^2 but not in C^3
"""
import numpy as np
import matplotlib.pyplot as plt

xmin, xmax = -1.0, 1.0

def fun(x):
    return x + np.exp(-np.abs(x)**3)
    
def dfun(x):
    return 1.0 - 3.0 * np.exp(-np.abs(x)**3) * x**2 * np.sign(x)

npt = [51,101,201,401,801]
errb, errf, errc, hh = [], [], [], []
for n in npt:
    x = np.linspace(xmin,xmax,n)
    h = (xmax - xmin)/(n-1); hh.append(h)
    f = fun(x)
    # backward difference
    dfb = (f[1:-1] - f[0:-2])/h
    # forward difference
    dff = (f[2:] - f[1:-1])/h
    # central difference
    dfc = (f[2:] - f[0:-2])/(2*h)
    # exact derivative
    dfe= dfun(x[1:-1])
    # errors
    plt.semilogy(x[1:-1], np.abs(dfe-dfc), label=str(n))
    errb.append(np.abs(dfe - dfb).max())
    errf.append(np.abs(dfe - dff).max())
    errc.append(np.abs(dfe - dfc).max())
    
plt.xlabel('x')
plt.ylabel('Error of central difference scheme')
plt.legend()

# Compute convergence rates
print("%5s %17s %10s %10s" % ('h','Backward','Forward','Central'))
for i in range(1,len(hh)):
    pb = np.log(errb[i-1]/errb[i])/np.log(2)
    pf = np.log(errf[i-1]/errf[i])/np.log(2)
    pc = np.log(errc[i-1]/errc[i])/np.log(2)
    print("%e   %f   %f   %f" % (hh[i],pb,pf,pc))

# Plot error vs h
plt.figure()
plt.loglog(hh,errb,'o-',hh,errf,'*-',hh,errc,'s--')
plt.xlabel('h')
plt.ylabel('Maximum error')
plt.legend(('Backward','Forward','Central'),loc='lower right')
plt.title('First derivative for $u(x)=x+\exp(-|x|^3)$, $x\in[-1,1]$')
plt.show()
