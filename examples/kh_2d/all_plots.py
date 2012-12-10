#!/usr/bin/env python

import glob, subprocess, sys

def run_with_output(argv):
  print 'running: ', argv
  subprocess.call(argv, stdout = sys.stdout, stderr = sys.stderr)

if __name__ == '__main__':
  run_with_output('./plot_energy.py')

  movieprops = [
    ('ne', 'density', 'Reds'),
    ('ni', 'density', 'Blues'),
    ('potential', 'potential', 'RdBu_r'),
    ('ex', 'potential', 'RdBu_r'),
    ('ey', 'potential', 'RdBu_r'),
    ('vex', 'potential', 'RdBu_r'),
    ('vey', 'potential', 'RdBu_r'),
    ('vix', 'potential', 'RdBu_r'),
    ('viy', 'potential', 'RdBu_r')
  ]

  for (fieldname, fieldtype, cmap) in movieprops:
    argv = ['./movie_heatmap.py', fieldtype, cmap]
    filenames = glob.glob('fields/{0}_??????.bin'.format(fieldname))
    argv.extend(filenames)
    run_with_output(argv)
    
    argv = ['./movie_yaverage.py', fieldtype]
    argv.extend(filenames)
    run_with_output(argv)
