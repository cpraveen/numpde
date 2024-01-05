'''
Generates mesh exterior to ellipse by conformal map
'''
import numpy as np
import matplotlib.pyplot as plt

# Ellipse (x/a)^2 + (y/b)^2 = 1
# ratio = a/b
ratio = 6.0
b = 1.0/np.sqrt(ratio**2 - 1.0)
a = ratio * b

print("Axis a = ", a)
print("Axis b = ", b)

n = 100 # no. angular points

R1, R2 = a+b, 100.0

m = (n/(2*np.pi)) * np.log(R2/R1)
m = int(m)

# Map uniform partition of xi in [0,1] to r in [R1,R2]
xi = np.linspace(0.0,1.0,m)
r = R1 * np.exp(xi * np.log(R2/R1))
print("Points in angular direction = ", n)
print("Points in radial  direction = ", m)

# Angular grid
theta = np.linspace(0.0,2.0*np.pi,n)

# Grid exterior to circle
# Cartesian coordinates
x = np.outer(r, np.cos(theta))
y = np.outer(r, np.sin(theta))

# Map circles to ellipses
z = x + 1j * y
z = 0.5 * ( z + 1/z)
x = np.real(z)
y = np.imag(z)

# inner boundary = ellipse
xe, ye = x[0,:], y[0,:]
ratio = np.max(xe)/np.max(ye)

# outer boundary
xo, yo = x[-1,:], y[-1,:]
rout = np.min( np.sqrt(xo**2 + yo**2) )
print("Distance of outer boundary = ", rout)
print("rout/a = ", rout/a)

plt.figure()
plt.plot(xe, ye)
plt.title('Ellipse, aspect ratio = '+str(ratio))
plt.axis('equal')

plt.figure()
# Plot angular lines
for i in range(m):
    plt.plot(x[i,:],y[i,:],'k-')
# Plot radial lines
for j in range(n):
    plt.plot(x[:,j],y[:,j],'k-')
plt.title(str(n)+' x '+str(m)+' mesh')
plt.axis('equal')
plt.show()
