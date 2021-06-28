# 1D Euler using FVM

* euler_fo: First order method, Lax-Friedrich flux
* euler_ho_1: Higher order, Lax-Friedrich flux
* euler_ho_2: Higher order, many fluxes

Compile using the makefile

```shell
make
```

which compiles many executables

## euler_fo

```shell
./euler_fo
gnuplot sod.gnu
open sol.pdf
```

## euler_ho_1

Set some parameters in `input.txt` file.

```shell
./euler_ho
gnuplot sod.gnu
open sol.pdf
```

## euler_ho_2

Of the three codes, this has the most variety of schemes. See the `input.txt` file for available options. The test cases are defined in files `sod.F90` and `shuosher.F90`. You can add a new test case and include the file.


### Sod test case

Set some parameters in `input.txt` file.

```shell
./euler_sod
gnuplot sod.gnu
open sol.pdf
```

### Shu-Osher test case

Set some parameters in `input.txt` file.

```shell
./euler_shuosher
gnuplot shuosher.gnu
open sol.pdf
```

Compare the results from different fluxes, reconstruction schemes and variable used for reconstruction.