push!(LOAD_PATH,".")
using VTE2d 
using Printf
using PyPlot

# Domain is unit square
xmin, xmax = 0.0, 1.0
ymin, ymax = xmin, xmax

# dx = dy
nx = ny = 20
h = (xmax-xmin)/(nx-1)

# Mesh
x = LinRange(xmin,xmax,nx)
y = LinRange(ymin,ymax,ny)
X, Y = meshgrid(x,y)

# Solution variables
psi   = zeros(nx,ny)
omega = zeros(nx,ny)

omega[:,:] = @. sin(2 * π * X) * sin(2 * π * Y)

res0, res, iter = solve_stream(psi, omega, h)
@printf("res0, res, iter = %e %e %d\n", res0, res, iter)
figure(figsize=(5,5))
contour(x, y, psi, 30)
show()
