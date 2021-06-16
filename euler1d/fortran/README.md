# 1D Euler using FVM

* euler_fo: First order method, Lax-Friedrich flux
* euler_ho_1: Higher order, Lax-Friedrich flux
* euler_ho_2: Higher order, many fluxes

Compile using the makefile

```bash
make
```

## euler_ho_2

Of the three codes, this has the more schemes. See the main function for available options.

Set the parameters in the main function, save, compile and run it

```bash
make euler_ho_2
./euler_ho_2
gnuplot shuosher.gnu
open sol.pdf
```

Compare the results from different fluxes, reconstruction schemes and variable used for reconstruction.