# 2-D BVP

## Poisson eqn in unit square: poisson1.cc

Compile the code
```
$ make poisson1
```
This creates the executable ```poisson1```. Run it
```
$ ./poisson1
```
Plot the solution in gnuplot as a surface plot
```
$ gnuplot
gnuplot> load 'splot.gnu'
```
Use your mouse to rotate the figure. You can also open the ```splot.eps``` file. Or plot filled colour plot
```
gnuplot> load 'contour.gnu'
```
and see the file ```contour.eps```.

## Same example using Eigen: poisson2.cc
