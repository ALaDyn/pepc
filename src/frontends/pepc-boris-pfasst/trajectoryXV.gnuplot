#!/usr/bin/gnuplot --persist

set term wx size 1000,1000
set output

set multiplot layout 2,3 rowsfirst

set label 1 'X' at graph 0.1,0.9 font ',14'
plot "fort.47" using 1:3 with l  title 'Boris classical',   \
     "fort.48" using 1:3 with l  title 'Boris SDC',         \
     "fort.49" using 1:3 with l  title 'Analytic solution', \
     "fort.50" using 1:3 with l  title 'Cyclotronic',       \
     "fort.51" using 1:3 with l  title 'Boris Patacchini',  \
     "fort.52" using 1:3 with l  title 'Boris leap-frog',   \
     "fort.53" using 1:3 with l  title 'Boris tan(alpha)/alpha',    \
     "fort.54" using 1:3 with l  title 'Tajima leap-frog implicit', \
     "fort.55" using 1:3 with l  title 'Tajima leap-frog explicit', \
     "fort.56" using 1:3 with l  title 'Inverting Verlet'
set label 1 'Y' at graph 0.1,0.9 font ',14'
plot "fort.47" using 1:4 with l  title 'Boris classical',   \
     "fort.48" using 1:4 with l  title 'Boris SDC',         \
     "fort.49" using 1:4 with l  title 'Analytic solution', \
     "fort.50" using 1:4 with l  title 'Cyclotronic',       \
     "fort.51" using 1:4 with l  title 'Boris Patacchini',  \
     "fort.52" using 1:4 with l  title 'Boris leap-frog',   \
     "fort.53" using 1:4 with l  title 'Boris tan(alpha)/alpha',    \
     "fort.54" using 1:4 with l  title 'Tajima leap-frog implicit', \
     "fort.55" using 1:4 with l  title 'Tajima leap-frog explicit', \
     "fort.56" using 1:4 with l  title 'Inverting Verlet'
set label 1 'Z' at graph 0.1,0.9 font ',14'
plot "fort.47" using 1:5 with l  title 'Boris classical',   \
     "fort.48" using 1:5 with l  title 'Boris SDC',         \
     "fort.49" using 1:5 with l  title 'Analytic solution', \
     "fort.50" using 1:5 with l  title 'Cyclotronic',       \
     "fort.51" using 1:5 with l  title 'Boris Patacchini',  \
     "fort.52" using 1:5 with l  title 'Boris leap-frog',   \
     "fort.53" using 1:5 with l  title 'Boris tan(alpha)/alpha',    \
     "fort.54" using 1:5 with l  title 'Tajima leap-frog implicit', \
     "fort.55" using 1:5 with l  title 'Tajima leap-frog explicit', \
     "fort.56" using 1:5 with l  title 'Inverting Verlet'

set label 1 'VX' at graph 0.1,0.9 font ',14'
plot "fort.47" using 1:6 with l  title 'Boris classical',   \
     "fort.48" using 1:6 with l  title 'Boris SDC',         \
     "fort.49" using 1:6 with l  title 'Analytic solution', \
     "fort.50" using 1:6 with l  title 'Cyclotronic',       \
     "fort.51" using 1:6 with l  title 'Boris Patacchini',  \
     "fort.52" using 1:6 with l  title 'Boris leap-frog',   \
     "fort.53" using 1:6 with l  title 'Boris tan(alpha)/alpha',    \
     "fort.54" using 1:6 with l  title 'Tajima leap-frog implicit', \
     "fort.55" using 1:6 with l  title 'Tajima leap-frog explicit', \
     "fort.56" using 1:6 with l  title 'Inverting Verlet'
set label 1 'VY' at graph 0.1,0.9 font ',14'
plot "fort.47" using 1:7 with l  title 'Boris classical',   \
     "fort.48" using 1:7 with l  title 'Boris SDC',         \
     "fort.49" using 1:7 with l  title 'Analytic solution', \
     "fort.50" using 1:7 with l  title 'Cyclotronic',       \
     "fort.51" using 1:7 with l  title 'Boris Patacchini',  \
     "fort.52" using 1:7 with l  title 'Boris leap-frog',   \
     "fort.53" using 1:7 with l  title 'Boris tan(alpha)/alpha',    \
     "fort.54" using 1:7 with l  title 'Tajima leap-frog implicit', \
     "fort.55" using 1:7 with l  title 'Tajima leap-frog explicit', \
     "fort.56" using 1:7 with l  title 'Inverting Verlet'
set label 1 'VZ' at graph 0.1,0.9 font ',14'
plot "fort.47" using 1:8 with l  title 'Boris classical',   \
     "fort.48" using 1:8 with l  title 'Boris SDC',         \
     "fort.49" using 1:8 with l  title 'Analytic solution', \
     "fort.50" using 1:8 with l  title 'Cyclotronic',       \
     "fort.51" using 1:8 with l  title 'Boris Patacchini',  \
     "fort.52" using 1:8 with l  title 'Boris leap-frog',   \
     "fort.53" using 1:8 with l  title 'Boris tan(alpha)/alpha',    \
     "fort.54" using 1:8 with l  title 'Tajima leap-frog implicit', \
     "fort.55" using 1:8 with l  title 'Tajima leap-frog explicit', \
     "fort.56" using 1:8 with l  title 'Inverting Verlet'

unset multiplot
