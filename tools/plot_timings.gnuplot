#!/usr/bin/env gnuplot
#
# For this to work, it needs to be called from PEPC's root directory (the one
# containing 'src' and have timings in 'timings'.
# There should also be a frontend built that possibly contains user_timers

# If you want to compare two timers (on different axes), set the following
# to non-zero (the timer's id)
comparison_timer = 0

set term x11

set xlabel 'time/iterations'
set ylabel 'time / a.u.'

if (comparison_timer > 0) {
   set y2label 'time / a.u.'
   set ytics nomirror
   set y2tics nomirror
}
# Have array to store timer names
array timers[90]
# Create up-to-date list of timers
system "test -f /tmp/pepc.timers.gplt && rm /tmp/pepc.timers.gplt"
system "for id in \`seq 1 90\` ; do echo \"timers[$id]='undefined'\" >> /tmp/pepc.timers.gplt; done"
system "grep 'integer, parameter :: t_.*=.*' src/treecode/lpepcsrc/module_timings.f90 | awk '{print \"timers[\"$6\"]=\\\"\"$4\"\\\"\"}' >> /tmp/pepc.timers.gplt"
system "for fe in bin/* ; do fe_short=\`basename $fe\`; grep 'integer, parameter :: t_.*=.*' src/frontends/${fe_short}/module_user_timings.f90 | awk '{print \"timers[\"60+$8\"]=\\\"\"$4\"\\\"\"}' >> /tmp/pepc.timers.gplt; done"
load '/tmp/pepc.timers.gplt'

do for [col = 1:90] {
   colpt = col + 2
   compt = comparison_timer + 2

   stat [1e-20:] 'timing/timing_avg.dat' u colpt nooutput
   if (STATS_max > STATS_min) {
      print 'column ',col, timers[col]
      set title 'column '.col.': '.timers[col] noenhanced
      if (comparison_timer > 0) {
         plot 'timing/timing_avg.dat' u colpt w lp lc 1 t timers[col] noenhanced, \
              '' u compt axis x1y2 t timers[comparison_timer] noenhanced
      } else {
         plot 'timing/timing_avg.dat' u colpt w lp lc 1 t timers[col] noenhanced
      }
      pause -1
   }
}

# Select a few timers out of all
#do for [col in "8 11 12 16 18 20 41 42 44 61 62 63"] {
#   colpt = col + 2
#   lab = colpt - 2
#   print 'column ',col, timers[lab]
#   set title 'column '.col.': '.timers[lab] noenhanced
#   plot 'timing/timing_avg.dat' u colpt w lp lt 1 lc 1
#   pause -1
#}
system "rm /tmp/pepc.timers.gplt"
