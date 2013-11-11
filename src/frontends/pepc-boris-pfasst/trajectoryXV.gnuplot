#!/usr/bin/gnuplot --persist

set term wx size 1000,1000
set output

plot "fort.47" using 3:4 with lp title 'Boris classical', \
     "fort.48" using 3:4 with lp title 'Boris SDC'


set multiplot layout 2,3 rowsfirst

set label 1 'X' at graph 0.92,0.9 font ',8'
plot "fort.47" using 1:3 with l  title 'Boris classical',   \
     "fort.48" using 1:3 with l  title 'Boris SDC',         \
     "fort.49" using 1:3 with l  title 'Analytic solution', \
     "fort.50" using 1:3 with l  title 'Cyclotronic',       \
     "fort.51" using 1:3 with l  title 'Boris no tan-transformation'
set label 1 'Y' at graph 0.92,0.9 font ',8'
plot "fort.47" using 1:4 with l  title 'Boris classical',   \
     "fort.48" using 1:4 with l  title 'Boris SDC',         \
     "fort.49" using 1:4 with l  title 'Analytic solution', \
     "fort.50" using 1:4 with l  title 'Cyclotronic',       \
     "fort.51" using 1:4 with l  title 'Boris no tan-transformation'
set label 1 'Z' at graph 0.92,0.9 font ',8'
plot "fort.47" using 1:5 with l  title 'Boris classical',   \
     "fort.48" using 1:5 with l  title 'Boris SDC',         \
     "fort.49" using 1:5 with l  title 'Analytic solution', \
     "fort.50" using 1:5 with l  title 'Cyclotronic',       \
     "fort.51" using 1:5 with l  title 'Boris no tan-transformation'

set label 1 'VX' at graph 0.92,0.9 font ',8'
plot "fort.47" using 1:6 with l  title 'Boris classical',   \
     "fort.48" using 1:6 with l  title 'Boris SDC',         \
     "fort.49" using 1:6 with l  title 'Analytic solution', \
     "fort.50" using 1:6 with l  title 'Cyclotronic',       \
     "fort.51" using 1:6 with l  title 'Boris no tan-transformation'
set label 1 'VY' at graph 0.92,0.9 font ',8'
plot "fort.47" using 1:7 with l  title 'Boris classical',   \
     "fort.48" using 1:7 with l  title 'Boris SDC',         \
     "fort.49" using 1:7 with l  title 'Analytic solution', \
     "fort.50" using 1:7 with l  title 'Cyclotronic',       \
     "fort.51" using 1:7 with l  title 'Boris no tan-transformation'
set label 1 'VZ' at graph 0.92,0.9 font ',8'
plot "fort.47" using 1:8 with l  title 'Boris classical',   \
     "fort.48" using 1:8 with l  title 'Boris SDC',         \
     "fort.49" using 1:8 with l  title 'Analytic solution', \
     "fort.50" using 1:8 with l  title 'Cyclotronic',       \
     "fort.51" using 1:8 with l  title 'Boris no tan-transformation'

unset multiplot
