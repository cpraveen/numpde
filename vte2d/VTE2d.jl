module VTE2d

using Printf
using PyPlot

#------------------------------------------------------------------------------
# Tensor product grid
# X[i,j] = x[i], Y[i,j] = y[j]
#------------------------------------------------------------------------------
function meshgrid(x,y)
   nx = length(x)
   ny = length(y)
   X = Array(x) .* Array(ones(ny))'
   Y = Array(ones(nx)) .* Array(y)'
   return X, Y
end

#------------------------------------------------------------------------------
# L2 residual norm for stream function equation
#------------------------------------------------------------------------------
function residual_stream(psi, omega, h)

   nx, ny = size(psi)

   res = 0.0
   for i=2:nx-1
      for j=2:ny-1
         tmp = (  (psi[i-1,j] - 2.0 * psi[i,j] + psi[i+1,j])/h^2 
                + (psi[i,j-1] - 2.0 * psi[i,j] + psi[i,j+1])/h^2
                + omega[i,j])
         res = res + tmp^2
      end
   end

   return sqrt(res/((nx-2)*(ny-2)))

end

#------------------------------------------------------------------------------
# Solve for stream function
#------------------------------------------------------------------------------
function solve_stream!(psi, omega, h, RTOL=1.0e-6, ITMAX=1000)

   nx, ny = size(psi)
   r = 2.0/(1.0 + π * h) # SOR factor

   res0 = residual_stream(psi, omega, h)
   res  = res0
   iter = 0

   while res > RTOL * res0 && iter < ITMAX
      # Update only interior values
      for i=2:nx-1
         for j=2:ny-1
            tmp = 0.25*(psi[i-1,j] + psi[i+1,j] + psi[i,j-1] + psi[i,j+1]
                        + h^2 * omega[i,j])
            psi[i,j] = (1.0 - r) * psi[i,j] + r * tmp
         end
      end

      res = residual_stream(psi, omega, h)
      iter += 1
   end

   if res > RTOL * res0 && iter == ITMAX
      println("Error: no convergence in stream function")
      println("Error: increase RTOL and/or ITMAX")
      @printf("res0, res, iter = %e %e %d\n", res0, res, iter)
      exit()
   end

   return res0, res, iter
end

#------------------------------------------------------------------------------
# Compute velocity at interior points
# Boundary point velocity are assumed not to change with time
#------------------------------------------------------------------------------
function compute_velocity!(u, v, psi, h)
   nx, ny = size(u)

   for i=2:nx-1
      for j=2:ny-1
         u[i,j] =  (psi[i,j+1] - psi[i,j-1]) / (2*h)
         v[i,j] = -(psi[i+1,j] - psi[i-1,j]) / (2*h)
      end
   end
end

#------------------------------------------------------------------------------
# Update vorticity to new time level
#------------------------------------------------------------------------------
function update_vort!(omega, psi, u, v, h, dt, nu)
   nx, ny = size(omega)
   w = copy(omega)

   # Interior points
   for i=2:nx-1
      for j=2:ny-1
         wx = (w[i+1,j] - w[i-1,j]) / (2*h)
         wy = (w[i,j+1] - w[i,j-1]) / (2*h)
         wxx = (w[i-1,j] - 2 * w[i,j] + w[i+1,j]) / h^2
         wyy = (w[i,j-1] - 2 * w[i,j] + w[i,j+1]) / h^2
         omega[i,j] = w[i,j] + dt*( - u[i,j] * wx - v[i,j] * wy + nu * (wxx + wyy) )
      end
   end

   # Bottom boundary
   for i=1:nx
      u1 = u[i,1]
      u2 = (psi[i,2] - psi[i,1]) / h
      u3 = (psi[i,3] - psi[i,2]) / h
      omega[i,1] = -(- 8 * u1 + 9 * u2 - u3) / (3*h)
   end

   # Top boundary
   for i=1:nx
      u1 = u[i,ny]
      u2 = (psi[i,ny] - psi[i,ny-1]) / h
      u3 = (psi[i,ny-1] - psi[i,ny-2]) / h
      omega[i,ny] = -( 8 * u1 - 9 * u2 + u3) / (3*h)
   end

   # Left boundary
   for j=2:ny-1
      v1 = v[1,j]
      v2 = -(psi[2,j] - psi[1,j]) / h
      v3 = -(psi[3,j] - psi[2,j]) / h
      omega[1,j] = (- 8 * v1 + 9 * v2 - v3) / (3*h)
   end

   # Right boundary
   for j=2:ny-1
      v1 = v[nx,j]
      v2 = -(psi[nx,j] - psi[nx-1,j]) / h
      v3 = -(psi[nx-1,j] - psi[nx-2,j]) / h
      omega[nx,j] = ( 8 * v1 - 9 * v2 + v3) / (3*h)
   end
end

#------------------------------------------------------------------------------
# Main function which solves lid-driven cavity problem
#------------------------------------------------------------------------------
function solve(Re, n)

   nu = 1.0/Re

   Tf = 100.0
   ITMAX = 100

   # Domain is unit square
   xmin, xmax = 0.0, 1.0
   ymin, ymax = xmin, xmax

   # dx = dy
   nx = ny = n
   h = (xmax-xmin)/(nx-1)

   # Mesh
   x = LinRange(xmin,xmax,nx)
   y = LinRange(ymin,ymax,ny)

   # Meshgrid
   X, Y = meshgrid(x,y)

   # Solution variables
   psi   = zeros(nx,ny)
   omega = zeros(nx,ny)
   u     = zeros(nx,ny)
   v     = zeros(nx,ny)

   # Velocity at top lid
   u[:,ny] = @. sin(π * x)^2

   dt = 0.25 * Re * h^2

   t, iter = 0.0, 0
   figure(figsize=(6,6))
   while t < Tf && iter < ITMAX
      res0, res, it = solve_stream!(psi, omega, h)
      @printf("stream: res0,res,it    = %12.4e %12.4e %4d\n",res0,res,it)
      compute_velocity!(u, v, psi, h)
      update_vort!(omega, psi, u, v, h, dt, nu)
      t += dt; iter += 1
      @printf("vort  : iter,t,min,max = %4d %12.4e %12.4e %12.4e\n",
              iter,t,minimum(omega),maximum(omega))
      figure(1); clf()
      contour(X,Y,omega,20)
      title("Vorticity")
      draw(); pause(0.1)
   end
   compute_velocity!(u, v, psi, h)

   return X, Y, psi, omega, u, v
end

#------------------------------------------------------------------------------
export meshgrid
export solve_stream!
export solve
#------------------------------------------------------------------------------

end
