#!/usr/bin/gnuplot --persist

set term wx size 1000,1000
set output

plot "fort.47" using 3:4 with lp title 'Boris classical',   \
     "fort.48" using 3:4 with lp title 'Boris SDC',         \
     "fort.49" using 3:4 with l  title 'Analytic solution', \
     "fort.50" using 3:4 with lp title 'Cyclotronic',       \
     "fort.51" using 3:4 with l  title 'Boris Patacchini',  \
     "fort.52" using 3:4 with l  title 'Boris leap-frog' ,  \
     "fort.53" using 3:4 with l  title 'Boris tan(alpha)/alpha',    \
     "fort.54" using 3:4 with lp title 'Tajima leap-frog implicit', \
     "fort.55" using 3:4 with lp title 'Tajima leap-frog explicit', \
     "fort.56" using 3:4 with lp title 'Inverting Verlet', \
     "fort.57" using 3:4 with lp title 'Cyclotronic (no tan trans)', \
     "fort.58" using 3:4 with lp title 'Boris Patacchini (no tan trans)'
