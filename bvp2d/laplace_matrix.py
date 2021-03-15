from numpy import ones
from scipy.sparse import spdiags,eye,kron
nx, ny = 6, 6       # Mesh size
dx, dy = 1.0, 1.0
mx, my = nx-2, ny-2 # No of interior points
Dxx = (1.0/dx**2) * spdiags([ones(mx),-2.0*ones(mx),ones(mx)],[-1,0,1],mx,mx)
Dyy = (1.0/dy**2) * spdiags([ones(my),-2.0*ones(my),ones(my)],[-1,0,1],my,my)
Ix, Iy = eye(mx), eye(my)
A = - kron(Iy, Dxx) - kron(Dyy, Ix)

print("Sparse Dx ="); print(Dx)
print("Full   Dx ="); print(Dx.todense())
print("Sparse Ix ="); print(Ix)
print("Full   Ix ="); print(Ix.todense())
print("A         = "); print(A.todense())
