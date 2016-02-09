set xlabel 'x'
set ylabel 'y'
splot 'u.dat' w l
set term postscript enhanced color
set out 'splot.eps'
replot
set out
set term x11
