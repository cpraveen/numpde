# 1-D FV code for Euler equations

Available options

* Test cases: sod
* Numerical fluxes: KFVS, Roe
* Reconstruction: first, minmod/mc, vanleer

Compile the code

```shell
make
```

Run first order code

```shell
./fv sod kfvs
gnuplot fvm.gnu
open fvm.ps
```

Run high order code: minmod limiter

```shell
./fv sod kfvs minmod
gnuplot fvm.gnu
open fvm.ps
```

Run high order code: vanleer limiter

```shell
./fv sod kfvs vanleer
gnuplot fvm.gnu
open fvm.ps
```
