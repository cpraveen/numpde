set term pdf
set out 'sol.pdf'

set xlabel 'x'
set ylabel 'Density'
p 'sol.txt' u 1:2 w lp

set xlabel 'x'
set ylabel 'Velocity'
p 'sol.txt' u 1:3 w lp

set xlabel 'x'
set ylabel 'Pressure'
p 'sol.txt' u 1:4 w lp
