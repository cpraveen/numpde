import PyPlot; plt = PyPlot
using LaTeXStrings

n = 20
xmin, xmax = 0.0, 1.0

f(x) = 2.0 + sin(15*x) * exp(-x)

x = LinRange(xmin, xmax, 10*n)
xg = LinRange(xmin, xmax, n)

plt.figure()
plt.plot(x, f.(x), "-")
plt.ylim(-0.1,3)
plt.xlabel("x"); plt.ylabel("f(x)"); plt.grid(true)
plt.title("Some function f(x)")
plt.show()

plt.figure()
plt.plot(x, f.(x), "-")
plt.plot(xg, 0*xg, "o")
plt.ylim(-0.10,3)
plt.xlabel("x"); plt.ylabel("f(x)"); plt.grid(true)
plt.title("Make a grid of n points")
plt.show()

plt.figure()
for i in 1:n
    plt.plot([xg[i],xg[i]],[0,f(xg[i])],"k-")
end
plt.plot(x, f.(x), "-")
plt.plot(xg, f.(xg), "s")
plt.plot(xg, 0*xg, "o")
plt.ylim(-0.10,3)
plt.xlabel("x"); plt.ylabel("f(x)"); plt.grid(true)
plt.title("Evaluate function on the grid")
plt.show()

plt.figure()
plt.plot(xg, f.(xg), "s")
plt.plot(xg, 0*xg, "o")
plt.ylim(-0.10,3)
plt.xlabel("x"); plt.ylabel("f(x)"); plt.grid(true)
plt.title(L"Approximate f(x) by data on grid: $(x_i, f(x_i))$, $i=1,2,\ldots,n$") 
plt.show()
