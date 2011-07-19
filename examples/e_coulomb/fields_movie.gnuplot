#!/usr/bin/gnuplot

# for greek letters - unfortunately only works with ps-output
# set terminal postscript enhanced

print "Rudimentary Gnuplot script for visualizing"
print "1D representation of spherically symmetric"
print "field and particle data from flat ascii"
print "files with the name 'radial_fields.xxxxxx',"
print "where xxxxxx represents a six-digit-number,"
print "starting from 000000."
print "The script simply loops over all files and"
print "draws several plots."


# columns in radial_fields outputfiles:
# /-----------------------------------\
# |1|2  |3    |4      |5      |6  |7  |
# |r|n_i|rho_i|v_i_min|v_i_max|E_r|phi|
# \-----------------------------------/

# construct timestep counter
if (!exists("timestep")) timestep = 0

print sprintf("=========== Timestep %6.6d ==========", timestep)

# filename constructor
fieldsname(itime) = sprintf("radial_fields.%6.6d", itime)

# define function to generate rgb color value form input ranges 0.0..1.0
rgb(r,g,b) = int(255*r)*65536 + int(255*g)*256 + int(255*b)

# color mapping function for charges: negative: blue, positive: red, neutral: grey
chargecolor(q) = q<0 ? rgb(0.,0.,1.) : ( q>0 ? rgb(1.,0.,0.) : rgb(.3,.3,.3) )

# set aspect ratio to 1:1:1
set size 1,1

actfieldsname = fieldsname(timestep)

set multiplot
set size square 0.45, 0.45

set origin 0.025,0.0
set xrange[0:8]
set yrange[1e-5:1]
set xlabel "r"
set ylabel "n_i"
set logscale y
set style fill pattern 2 border 1
plot actfieldsname using 1:2:(0.) notitle with filledcurves above, '' u 1:2 with lines lt -1 title 'n_i' 

set origin 0.0,0.525
set xrange[0:8]
set yrange[1e-5:1]
set xlabel "r"
set ylabel "rho_i"
set logscale y
set style fill pattern 1 border 1
plot actfieldsname using 1:3:(0.) notitle with filledcurves above, '' u 1:3 with lines lt -1 title 'rho_i' 

set origin 0.525,0.0
set xrange[0:8]
set yrange[1e-5:15]
set xlabel "r"
set ylabel "v_i"
set logscale y
set style fill pattern 3 border 1
plot actfieldsname using 1:4:5 notitle with filledcurves, '' u 1:4 with lines lt -1 title 'v_i_min', '' u 1:5 with lines lt -1 title 'v_i_max'

set origin 0.525,0.525
set xrange[0:8]
#set yrange[1e-5:15]
set xlabel "r"
set ylabel "E_r & phi"
set logscale y
set style fill pattern 4 border 1
plot actfieldsname using 1:6 with lines title 'E_r', '' using 1:7 with lines title 'phi' 



timestep = timestep + 1;

# check wether new file exists. if not: restart from the beginning
actfieldsname = fieldsname(timestep)
fileexists = system(sprintf('if [ -e "%s" ]; then echo "1"; else echo "0"; fi', actfieldsname))
if (fileexists == 0) timestep = 0

 pause 0.2
 reread
