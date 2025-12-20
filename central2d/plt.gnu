reset
set term pdf
set out 'res.pdf'

set grid

set xlabel 'x'
set ylabel 'v'
p 'line0.dat' u 1:4 w l,'line.dat' u 1:4 w p pt 6

set xlabel 'x'
set ylabel 'p'
p 'line0.dat' u 1:5 w l,'line.dat' u 1:5 w p pt 6

set term x11
set out
