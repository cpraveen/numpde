set term postscript enhanced color
set out 'fvm.ps'

set title 'Density'
p 'sol.dat' u 1:2 t 'FVM' w p pt 6, \
  'msod_exact.dat' u 1:2 t 'Exact' w l lw 2

set title 'Velocity'
p 'sol.dat' u 1:3 t 'FVM' w p pt 6, \
  'msod_exact.dat' u 1:3 t 'Exact' w l lw 2

set title 'Pressure'
p 'sol.dat' u 1:4 t 'FVM' w p pt 6, \
  'msod_exact.dat' u 1:4 t 'Exact' w l lw 2
