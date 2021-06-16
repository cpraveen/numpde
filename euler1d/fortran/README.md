# 1D Euler using FVM

* euler_fo: First order method, Lax-Friedrich flux
* euler_ho_1: Higher order, Lax-Friedrich flux
* euler_ho_2: Higher order, many fluxes

Compile using the makefile

```bash
make
```

## euler_ho_2

Of the three codes, this has the most variety of schemes. See the main function for available options. The test cases are defined in files `sod.F90` and `shuosher.F90`. You can add a new test case and include the file.

Set the parameters in the main function, save, compile and run it to solve Shu-Osher test case.

```bash
make euler_ho_2
./euler_ho_2
gnuplot shuosher.gnu
open sol.pdf
```

Compare the results from different fluxes, reconstruction schemes and variable used for reconstruction.