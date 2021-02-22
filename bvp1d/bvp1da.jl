# Solve -u'' = f using finite difference
import PyPlot; plt = PyPlot
using LinearAlgebra

# RHS function
f(x) = sin(x)

# exact solution
uexact(x) = sin(x)

# Domain
xmin, xmax = 0.0, 2.0*pi

# Grid of n points
n = 20
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
println("Max error = ", maximum(abs.(uexact.(x)-u)))

# Exact solution on fine mesh for plotting
xe = LinRange(xmin, xmax, 100); ue = uexact.(xe)

# Plot exact and numerical solution
plt.plot(xe,ue,x,u,"o")
plt.legend(("Exact solution","Numerical solution"))
plt.xlabel("x")
plt.ylabel("u")
plt.show()
