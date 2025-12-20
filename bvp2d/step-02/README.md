# Sparse matrix, iterative solvers and BVP

Compile

```shell
make
```

How to use initialize sparse matrices

* `sparse_test1`: using row, col and val arrays
* `sparse_test2`: by setting individual elements

## Solve 1D Poisson equation

Get help

```shell
./linsol_test_1d
```

Run

```shell
./linsol_test_1d 50 sor 1000
```

Plot solution in gnuplot

```gnuplot
gnuplot> load 'plot.gnu'
```

## Solve 2D Poisson equation

Get help

```shell
./linsol_test_2d
```

Run

```shell
./linsol_test_2d 50 50 jacobi 5000
```

Plot solution in gnuplot

```gnuplot
gnuplot> load 'contour.gnu'
gnuplot> load 'splot.gnu'
```

and also see the files `contour.eps` and `splot.eps`.
