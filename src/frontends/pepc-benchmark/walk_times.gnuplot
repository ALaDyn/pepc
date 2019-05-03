#!/usr/bin/env gnuplot
set term x11
set term pdf size 21cm,30cm
set out 'walk_times.pdf'

!grep "tree walk time" log > walk_times.dat
!awk '{print $6}' walk_times.dat | gsl-histogram 0 3 300 > walk_times.histo

set multiplot title 'walk times vs time step' layout 2,1

set xlabel 'time step'
set ylabel 'time [s]'
plot 'walk_times.dat' u :6 w l title 'walk time [s]'

set style data histogram
set boxwidth 2.7 relative
set bars fullwidth
set style fill solid 2 border lt -1
set xtics nomirror ("0" 0, "0.5" 49, "1.0" 99, "1.5" 149, "2.0" 199, "2.5" 249, "3.0" 299)
set xlabel 'walk time [s]'
set ylabel 'occurance'
plot 'walk_times.histo' u 3 notitle lc 2

unset multiplot

#pause -1
