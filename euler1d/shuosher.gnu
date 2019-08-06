set term pdf
set out 'sol.pdf'

set xlabel 'x'
set ylabel 'Density, rho'
p 'sol.txt' u 1:2 t 'Numerical' w l, 'shuosher.txt' u 1:2 t 'Exact' w l

set xlabel 'x'
set ylabel 'Momentum, rho * u'
p 'sol.txt' u 1:($2*$3) t 'Numerical' w l, 'shuosher.txt' u 1:3 t 'Exact' w l

set xlabel 'x'
set ylabel 'Total energy, E'
p 'sol.txt' u 1:(($4)/0.4+0.5*$2*$3*$3) t 'Numerical' w l, 'shuosher.txt' u 1:4 t 'Exact' w l
