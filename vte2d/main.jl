push!(LOAD_PATH,".")
using VTE2d
using PyPlot

# Reynolds number
Re = 10.0
# Number of grid points
n = 20

X,Y,psi,omega,u,v = solve(Re, n)

# Plot final solution contours
figure(figsize=(10,5))
subplot(121)
contour(X,Y,psi,20); draw()
title("Stream function")
subplot(122)
contour(X,Y,omega,20)
title("Vorticity")
savefig("stream_vort.svg")

# Plot velocity along center line
figure(figsize=(10,5))
m = Int(n/2)
subplot(121)
plot(u[m,:], Y[m,:], "o-")
xlabel("u"); ylabel("y"); title("u(0.5,y)")
grid(true)
subplot(122)
plot(X[:,m], v[:,m], "s-")
xlabel("x"); ylabel("v"); title("v(x,0.5)")
grid(true)
savefig("velocity1.svg")

# Plot velocity vectors
figure(figsize=(5,5))
plot([0,1,1,0,0],[0,0,1,1,0],"b-",lw=0.5) # Domain
quiver(X,Y,u,v)
title("Velocity vectors")
savefig("velocity2.svg")

show()

