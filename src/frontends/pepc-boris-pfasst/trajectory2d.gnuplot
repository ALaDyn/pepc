#!/usr/bin/gnuplot --persist

set term wx size 1000,1000
set output

plot "fort.47" using 3:4 with lp title 'Boris classical',   \
     "fort.48" using 3:4 with lp title 'Boris SDC',         \
     "fort.49" using 3:4 with l  title 'Analytic solution', \
     "fort.50" using 3:4 with lp title 'Cyclotronic', \
     "fort.51" using 3:4 with l  title 'Boris no tan-transformation', \
     "fort.52" using 3:4 with l  title 'Boris leap-frog'     

