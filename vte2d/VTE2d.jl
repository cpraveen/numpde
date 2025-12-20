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
# Solve for stream function using SOR
#------------------------------------------------------------------------------
function solve_stream!(psi, omega, h, RTOL=1.0e-4, ITMAX=1000)

   nx, ny = size(psi)
   r = 2.0/(1.0 + π * h) # SOR factor

   psi .= 0.0
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
# Time step from Fourier stability
#------------------------------------------------------------------------------
function time_step(nu, h, u, v)
   nx, ny = size(u)

   dt = 0.25 * h^2 / nu

   for i=1:nx
      for j=1:ny
         dt = min(dt, 2.0*nu/(u[i,j]^2 + v[i,j]^2 + 1.0e-12))
      end
   end

   return dt
end

#------------------------------------------------------------------------------
# Fill vorticity at boundary points
#------------------------------------------------------------------------------
function boundary_vorticity!(omega, psi, u, v, h)
   nx, ny = size(omega)
   ih2 = 2.0 / h^2

   for i=2:nx-1
      # Bottom boundary
      omega[i,1] = ih2 * (psi[i,1] - psi[i,2] + h * u[i,1])
      # Top boundary
      omega[i,ny] = ih2 * (psi[i,ny] - psi[i,ny-1] - h * u[i,ny])
   end

   for j=2:ny-1
      # Left boundary
      omega[1,j] = ih2 * (psi[1,j] - psi[2,j] - h * v[1,j])
      # Right boundary
      omega[nx,j] = ih2 * (psi[nx,j] - psi[nx-1,j] + h * v[nx,j])
   end
end

#------------------------------------------------------------------------------
# Update vorticity to new time level
#------------------------------------------------------------------------------
function update_vort!(omega, psi, u, v, h, dt, nu)
   nx, ny = size(omega)
   w = copy(omega)

   # Interior points
   res = 0.0
   for i=2:nx-1
      for j=2:ny-1
         wx = (w[i+1,j] - w[i-1,j]) / (2*h)
         wy = (w[i,j+1] - w[i,j-1]) / (2*h)
         wxx = (w[i-1,j] - 2 * w[i,j] + w[i+1,j]) / h^2
         wyy = (w[i,j-1] - 2 * w[i,j] + w[i,j+1]) / h^2
         omega[i,j] = w[i,j] + dt*( - u[i,j] * wx - v[i,j] * wy + nu * (wxx + wyy) )
         res += ((omega[i,j] - w[i,j])/dt)^2
      end
   end

   return sqrt(res/((nx-2)*(ny-2)))
end

#------------------------------------------------------------------------------
# Main function which solves lid-driven cavity problem
#------------------------------------------------------------------------------
function solve(Re, n; Tf=50.0, ITMAX=50000, TOL=1.0e-4, pgap=100)

   nu = 1.0/Re

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
   u[:,ny] .= 1.0

   t, iter, wres = 0.0, 0, 1.0e20
   figure(figsize=(6,6))
   while t < Tf && iter < ITMAX && wres > TOL
      res0, res, it = solve_stream!(psi, omega, h)
      @printf("stream: res0,res,it    = %12.4e %12.4e %4d\n",res0,res,it)
      compute_velocity!(u, v, psi, h)
      dt = time_step(nu, h, u, v)
      boundary_vorticity!(omega, psi, u, v, h)
      wres = update_vort!(omega, psi, u, v, h, dt, nu)
      t += dt; iter += 1
      @printf("vort  : iter,t,res,min,max = %4d %10.3e %10.3e %10.3e %10.3e\n",
              iter,t,wres,minimum(omega),maximum(omega))
      if iter%pgap == 0
         figure(1); clf()
         contour(X,Y,psi,20)
         title("Streamlines")
         draw(); pause(0.1)
      end
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
