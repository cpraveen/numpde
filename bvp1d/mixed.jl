# Solve mixed problem using finite difference
# u'' = f in (0,1)
# u'(0) = sigma, u(1) = beta
using Printf
using PyPlot
using LinearAlgebra

# Domain
xmin, xmax = 0.0, 1.0

# RHS function
f(x) = exp(x)

# Boundary condition
sigma, beta = 0.0, 0.0

# exact solution
uexact(x) = 1.0 - x - exp(1.0) + exp(x)

# Grid of n points
function solve(n,method)
   h = (xmax - xmin)/(n - 1)
   x = LinRange(xmin,xmax,n)

   # array for solution
   u = zeros(n)

   # BC for last point
   u[end] = beta

   b = h^2 * f.(x[1:end-1])
   if method == 1
      b[1]  = h * sigma
   else
      b[1]  = h * sigma + 0.5 * h^2 * f(x[1])
   end
   b[end] -= beta

   m = n - 1 # exclude boundary conditions
   dl = ones(m-1)
   dm = -2.0 * ones(m); dm[1] = -1.0
   du = ones(m-1)
   A = Tridiagonal(dl, dm, du)
   u[1:end-1] = A \ b
   err = maximum(abs.(uexact.(x)-u))
   return h, x, u, err
end

function compare(n)
   # Exact solution for plotting
   xe = LinRange(xmin, xmax, 100)
   ue = uexact.(xe)

   h1, x1, u1, err1 = solve(n,1)
   h2, x2, u2, err2 = solve(n,2)
   figure()
   plot(x1,u1,"o",label="Method=1")
   plot(x2,u2,"s",label="Method=2")
   plot(xe,ue,label="Exact")
   legend()
end

function convergence(method)
   # Compute error norm for different meshes
   h, err = zeros(0), zeros(0)
   for n in [20,40,80,160,320]
      h1, x1, u1, err1 = solve(n,method)
      push!(h, h1)
      push!(err, err1)
   end

   # Compute convergence rate in L2 norm
   println("------ Method = $method ------ ")
   @printf("h = %e   err = %e\n", h[1], err[1])
   for i in 2:length(h)
      p = log(err[i-1]/err[i])/log(2)
      @printf("h = %e   err = %e  rate = %f\n", h[i], err[i], p)
   end
   return h, err
end

compare(20)
h1, err1 = convergence(1)
h2, err2 = convergence(2)

figure()
loglog(h1, err1, "o-", label="Method=1")
loglog(h2, err2, "s-", label="Method=2")
xlabel("h")
ylabel("Error norm")
legend()
show()
