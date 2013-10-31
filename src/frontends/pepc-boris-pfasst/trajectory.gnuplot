#!/usr/bin/gnuplot --persist

set term aqua size 1000,1000
set output

splot "fort.47" using 3:4:5 with lines
