# 1-D FV code for Euler equations

Compile the code

```shell
make
```

Run first order code

```shell
./fv sod
gnuplot fvm.gnu
open fvm.ps
```

Run high order code: minmod limiter

```shell
./fv sod minmod
gnuplot fvm.gnu
open fvm.ps
```

Run high order code: vanleer limiter

```shell
./fv sod vanleer
gnuplot fvm.gnu
open fvm.ps
```