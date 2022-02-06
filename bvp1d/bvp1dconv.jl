# Solve -u'' = f using finite difference
using Printf
using PyPlot
using LinearAlgebra

# RHS function
f(x) = sin(x)

# exact solution
uexact(x) = sin(x)

# Domain
xmin, xmax = 0.0, 2.0*pi

# Grid of n points
function error(n)
   h = (xmax - xmin)/(n - 1)
   x = LinRange(xmin,xmax,n)

   # array for solution
   u = zeros(n)

   # BC for first and last points
   u[1]  = uexact(x[1])
   u[end] = uexact(x[end])

   b = h^2 * f.(x[2:end-1])
   b[1]  += u[1]
   b[end] += u[end]

   m = n - 2 # exclude boundary conditions
   A = Tridiagonal(-ones(m-1), 2*ones(m), -ones(m-1))
   u[2:end-1] = A \ b
   return h, maximum(abs.(uexact.(x)-u))
end

# Compute error norm for different meshes
h, err = zeros(0), zeros(0)
for n in [20,40,80,160,320,640,1280]
   h1, err1 = error(n)
   push!(h, h1)
   push!(err, err1)
end

# Compute convergence rate in L2 norm
@printf("h = %e   err = %e\n", h[1], err[1])
for i in 2:length(h)
   p = log(err[i-1]/err[i])/log(2)
   @printf("h = %e   err = %e  rate = %f\n", h[i], err[i], p)
end

loglog(h, err, "o-")
xlabel("h")
ylabel("Error norm")
show()
