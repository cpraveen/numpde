# 2-D BVP

## Poisson eqn in unit square: poisson1.cc

Compile the code
```
$ make poisson1
```
This creates the executable `poisson1`. Run it
```
$ ./poisson1
```
Plot the solution in gnuplot as a surface plot
```
$ gnuplot
gnuplot> load 'splot.gnu'
```
Use your mouse to rotate the figure. You can also open the `splot.eps` file. Or plot filled colour plot
```
gnuplot> load 'contour.gnu'
```
and see the file `contour.eps`.

To compile in debug which does some safety checks, compile like this
```
make poisson1 debug=yes
```

## Same example using Eigen: poisson2.cc

Use a file diff tool like opendiff (on macOS) or meld (on Linux) to see the differences between `poisson1.cc` and `poisson2.cc` files, or use diff in the terminal
```
diff poisson1.cc poisson2.cc
```

Compile and run in the same way as the previous code.

## Exercise

The code assumes that dx = dy, i.e., the grid point spacing is same in both directions. Modify it so that you can use different grid point spacing in the two directions.
