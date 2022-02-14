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
   ue = uexact.(x)
   du = ue - u
   err_max = maximum(abs.(du))
   # solution norm
   err_l2  = 0.5*h*du[1]^2 + h*sum(du[2:end-1].^2) + 0.5*h*du[end]^2
   err_l2  = sqrt(err_l2)
   # derivative norm
   deriv   = (du[2:end] - du[1:end-1])/h
   err_h1  = sqrt(h * sum(deriv.^2))
   return h, x, u, err_max, sqrt(err_l2^2 + err_h1^2)
end

function compare(n)
   # Exact solution for plotting
   xe = LinRange(xmin, xmax, 100)
   ue = uexact.(xe)

   h1, x1, u1, err1_max, err1_l2 = solve(n,1)
   h2, x2, u2, err2_max, err2_l2 = solve(n,2)
   figure()
   plot(x1,u1,"o",label="Method=1")
   plot(x2,u2,"s",label="Method=2")
   plot(xe,ue,label="Exact")
   xlabel("x"); ylabel("Error"); title(string(n," points"))
   legend()
end

function convergence(method)
   # Compute error norm for different meshes
   h, err_max, err_l2 = zeros(0), zeros(0), zeros(0)
   for n in [20,40,80,160,320]
      h1, x1, u1, err1_max, err1_l2 = solve(n,method)
      push!(h, h1)
      push!(err_max, err1_max)
      push!(err_l2,  err1_l2)
   end

   # Compute convergence rate in L2 norm
   println("------ Method = $method ------ ")
   #@printf("h = %e   err = %e\n", h[1], err_max[1])
   @printf("%20s max norm %20s Sobolev norm\n", " "," ")
   for i in 2:length(h)
      pmax = log(err_max[i-1]/err_max[i])/log(2)
      pl2  = log(err_l2[i-1]/err_l2[i])/log(2)
      @printf("h= %e  err= %e  rate= %4.2f  err= %e  rate= %4.2f\n", 
              h[i], err_max[i], pmax, err_l2[i], pl2)
   end
   return h, err_max, err_l2
end

compare(20)
h1, err1_max, err1_l2 = convergence(1)
h2, err2_max, err2_l2 = convergence(2)

figure()
loglog(h1, err1_max, "o-", label="Method=1")
loglog(h2, err2_max, "s-", label="Method=2")
xlabel("h")
ylabel("Max error norm")
legend()

figure()
loglog(h1, err1_l2, "o-", label="Method=1")
loglog(h2, err2_l2, "s-", label="Method=2")
xlabel("h")
ylabel("Sobolev error norm")
legend()

show()
