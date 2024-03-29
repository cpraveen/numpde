"""
Solve 2d poisson using sparse, direct solver
     -Laplace(u) = f
              u  = g
"""
from numpy import sin,cos,pi,linspace,meshgrid,zeros,ones,reshape,abs
from scipy.sparse import spdiags,eye,kron
from scipy.sparse.linalg import spsolve
import matplotlib.pyplot as plt

# Right hand side function
f = lambda x,y: sin(2.0*pi*x) * cos(2.0*pi*y)

# Exact solution
uexact = lambda x,y: 1.0/(2.0*(2.0*pi)**2) * sin(2.0*pi*x) * cos(2.0*pi*y)

# Domain
xmin, xmax = 0.0, 1.0
ymin, ymax = 0.0, 1.0

# Mesh size
nx, ny = 50, 50

dx, dy = (xmax-xmin)/(nx-1), (ymax-ymin)/(ny-1)
x = linspace(xmin,xmax,nx)
y = linspace(ymin,ymax,ny)
X,Y = meshgrid(x,y,indexing='ij')

print("Grid size = ", nx, "x", ny)
print("dx, dy    = ", dx, dy)

# No of interior points
mx, my = nx-2, ny-2

Dxx = (1.0/dx**2) * spdiags([ones(mx),-2.0*ones(mx),ones(mx)],[-1,0,1],mx,mx)
Dyy = (1.0/dy**2) * spdiags([ones(my),-2.0*ones(my),ones(my)],[-1,0,1],my,my)
Ix, Iy = eye(mx), eye(my)
A = - kron(Iy, Dxx) - kron(Dyy, Ix)

# Exact solution
uex = uexact(X,Y)

# RHS vector
b = f(X[1:-1,1:-1], Y[1:-1,1:-1])

# Add boundary values to RHS vector
b[ 0, :] += uex[   0,1:-1] / dx**2 # left   boundary
b[-1, :] += uex[  -1,1:-1] / dx**2 # right  boundary
b[ :, 0] += uex[1:-1,   0] / dy**2 # bottom boundary
b[ :,-1] += uex[1:-1,  -1] / dy**2 # top    boundary

b = reshape(b, mx*my, order='F')

# Solve
sol = spsolve(A,b)

# Reshape to array
u = uex.copy() # Fill boundary values
u[1:-1,1:-1] = reshape(sol,(mx,my),order='F')
print('Max error = ', abs(u-uex).max())

# Contour plot solution
plt.figure(figsize=(5,5))
plt.title("Solution")
plt.contour(X,Y,u,20)
plt.xlabel("x"); plt.ylabel("y")

# Color plot error
plt.figure()
plt.title("Error")
cs = plt.contourf(X,Y,abs(u-uex),levels=30)
plt.colorbar(cs)
plt.xlabel("x"); plt.ylabel("y")

plt.show()
