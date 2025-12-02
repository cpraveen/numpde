set xlabel 'x'
set ylabel 'u'
set grid
p 'burg2_ROE.txt' u 1:2 w lp lw 2 t 'ROE', \
  'burg2_LF.txt'  u 1:2 w lp lw 2 t 'LF', \
  'burg2_LLF.txt' u 1:2 w lp lw 2 t 'LLF', \
  'burg2_GLF.txt' u 1:2 w lp lw 2 t 'GLF', \
  'burg2_GOD.txt' u 1:2 w lp lw 2 t 'GOD', \
  'burg2_ROE.txt' u 1:3 w l  lw 2 t 'Exact'
