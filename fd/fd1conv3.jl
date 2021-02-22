using PyPlot

xmin, xmax = 0.0, 1.0

f(x)  = sin(2*pi*x)
df(x) = 2*pi*cos(2*pi*x)

dx, errb, errf, errc = zeros(0), zeros(0), zeros(0), zeros(0)
for n in [50,100,200,400,800]
   x = Array(LinRange(xmin,xmax,n)) # cannot modify LinRange, so make it Array
   h = (xmax - xmin)/(n-1)
   x[2:end-1] += 0.1 * h * (2*rand(n-2) .- 1)
   F = f.(x)
   i = 2:n-1 # Interior points
   dfb = (F[i] - F[i.-1])./(x[i] - x[i.-1])       # backward difference
   dff = (F[i.+1] - F[i])./(x[i.+1] - x[i])       # forward difference
   dfc = (F[i.+1] - F[i.-1])./(x[i.+1] - x[i.-1]) # central difference
   dfe = df.(x[i])
   push!(dx,maximum(x[2:end].-x[1:end-1]))
   push!(errb, maximum(abs.(dfb.-dfe)))
   push!(errf, maximum(abs.(dff.-dfe)))
   push!(errc, maximum(abs.(dfc.-dfe)))
end

loglog(dx, errb, "o-", label="Backward")
loglog(dx, errf, "*-", label="Forward")
loglog(dx, errc, "s--", label="Central")
xlabel("h"); ylabel("Error norm"); legend()
show()
