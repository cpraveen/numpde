# 1-D FV code for Euler equations

Compile the code

```shell
make
```

Run first order code

```shell
./fv 1
gnuplot fvm.gnu
open fvm.ps
```

Run high order code: minmod limiter

```shell
./fv 1 minmod
gnuplot fvm.gnu
open fvm.ps
```

Run high order code: vanleer limiter

```shell
./fv 1 vanleer
gnuplot fvm.gnu
open fvm.ps
```