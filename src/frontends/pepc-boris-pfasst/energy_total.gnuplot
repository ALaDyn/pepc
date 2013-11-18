#!/usr/bin/gnuplot --persist

set term wx size 1000,1000
set output

set multiplot layout 3,3 rowsfirst

plot "fort.147" using 1:4 with l  title 'Boris classical'
plot "fort.148" using 1:4 with l  title 'Boris SDC'
plot "fort.149" using 1:4 with l  title 'Analytic solution'
plot "fort.150" using 1:4 with l  title 'Cyclotronic'
plot "fort.151" using 1:4 with l  title 'Boris Patacchini'
plot "fort.152" using 1:4 with l  title 'Boris leap-frog'
plot "fort.153" using 1:4 with l  title 'Boris tan(alpha)/alpha'
plot "fort.154" using 1:4 with l  title 'Tajima leap-frog implicit'
plot "fort.155" using 1:4 with l  title 'Tajima leap-frog explicit'
