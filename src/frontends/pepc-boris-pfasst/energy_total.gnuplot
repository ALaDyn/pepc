#!/usr/bin/gnuplot --persist

set term wx size 1000,1000
set output

f(x) = mean_y

set multiplot layout 3,4 rowsfirst

fit f(x) "fort.147" u 1:4 via mean_y
plot "fort.147" using 1:($4/mean_y) with l  title 'Boris classical'

fit f(x) "fort.148" u 1:4 via mean_y
plot "fort.148" using 1:($4/mean_y) with l  title 'Boris SDC'

fit f(x) "fort.149" u 1:4 via mean_y
plot "fort.149" using 1:($4/mean_y) with l  title 'Analytic solution'

fit f(x) "fort.150" u 1:4 via mean_y
plot "fort.150" using 1:($4/mean_y) with l  title 'Cyclotronic'

fit f(x) "fort.151" u 1:4 via mean_y
plot "fort.151" using 1:($4/mean_y) with l  title 'Boris Patacchini'

fit f(x) "fort.152" u 1:4 via mean_y
plot "fort.152" using 1:($4/mean_y) with l  title 'Boris leap-frog'

fit f(x) "fort.153" u 1:4 via mean_y
plot "fort.153" using 1:($4/mean_y) with l  title 'Boris tan(alpha)/alpha'

fit f(x) "fort.154" u 1:4 via mean_y
plot "fort.154" using 1:($4/mean_y) with l  title 'Tajima leap-frog implicit'

fit f(x) "fort.155" u 1:4 via mean_y
plot "fort.155" using 1:($4/mean_y) with l  title 'Tajima leap-frog explicit'

fit f(x) "fort.156" u 1:4 via mean_y
plot "fort.156" using 1:($4/mean_y) with l  title 'Inverting Verlet'

fit f(x) "fort.157" u 1:4 via mean_y
plot "fort.157" using 1:($4/mean_y) with l  title 'Cyclotronic (no tan trans)'

fit f(x) "fort.158" u 1:4 via mean_y
plot "fort.158" using 1:($4/mean_y) with l  title 'Boris Patacchini (no tan trans)'

