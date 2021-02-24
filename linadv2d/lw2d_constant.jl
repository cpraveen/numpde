using ShiftedArrays
using PyPlot
using Printf

nx, ny = 100, 100      # grid size
xmin, xmax = -1.0, 1.0 # domain size
ymin, ymax = -1.0, 1.0
u, v = 1.0, 1.0        # advection velocity
Tf = 2.0               # one time period

# make grid
dx, dy = (xmax - xmin)/nx, (ymax - ymin)/ny
x = LinRange(xmin+0.5*dx, xmax-0.5*dx, nx); x = x'
y = LinRange(ymin+0.5*dy, ymax-0.5*dy, ny)

# set initial condition
q = @. sin(2*pi*x) * sin(2*pi*y)

# save copy of ic, we use it to compute error norm
q0 = q

# plot initial condition
contour(x,y,q)
xlabel("x"); ylabel("y"); axis("equal")
title("Initial condition")
draw(); pause(1)

# dt from cfl condition
dt = 0.72/(abs(u)/dx + abs(v)/dy)
println("Grid size nx, ny = ", nx, ny)
println("Grid size dx, dy = ", dx, dy)
println("Time step        = ", dt)

s1, s2 = u*dt/dx, v*dt/dy

# Time loop
t, it = 0.0, 0
while t < Tf
   global dt, q, t, it
   if t+dt > Tf # adjust dt so we reach Tf exactly
       dt = Tf - t
   end
   # generate stencil values using periodicity
   qim1j = CircShiftedArray(q,(1,0))
   qip1j = CircShiftedArray(q,(-1,0))
   qijm1 = CircShiftedArray(q,(0,1))
   qijp1 = CircShiftedArray(q,(0,-1))
   qim1jm1 = CircShiftedArray(q,(1,1))
   qim1jp1 = CircShiftedArray(q,(1,-1))
   qip1jm1 = CircShiftedArray(q,(-1,1))
   qip1jp1 = CircShiftedArray(q,(-1,-1))
   # LW scheme
   q = (q - 0.5 * s1 * (qip1j - qim1j) - 0.5 * s2 * (qijp1 - qijm1)
         + 0.5 * s1^2 * (qim1j - 2.0*q + qip1j)
         + 0.25 * s1 * s2 * (qip1jp1 - qip1jm1 - qim1jp1 + qim1jm1)
         + 0.5 * s2^2 * (qijm1 - 2.0*q + qijp1))
   t += dt; it += 1
   println("it,t = ",it,"  ",t)
   if mod(it,5) == 0 || abs(t-Tf) < 1.0e-12
       clf()
       contour(x,y,q)
       title(string("t = ",t))
       xlabel("x"); ylabel("y"); axis("equal")
       draw(); pause(0.1)
   end
end

# Compute error norm: initial condition is exact solution
l1_err = sum(abs.(q-q0)) / (nx*ny)
l2_err = sqrt(sum((q-q0).^2) / (nx*ny))
li_err = maximum(abs.(q-q0))
@printf("dx,dy,l1,l2,linf error = %10.4e %10.4e %10.4e %10.4e %10.4e\n",
        dx,dy,l1_err,l2_err,li_err)
show()
