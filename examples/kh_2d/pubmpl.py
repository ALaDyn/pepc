#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

def set_params(relwidth = 0.9, labelsize_pt = 10, ticksize_pt = 8, columnwidth_pt = 246, ppi = 72.27, aspect_ratio = (1.0 + np.sqrt(5)) / 2.0):
  columnwidth = columnwidth_pt / ppi
  fig_width = relwidth * columnwidth
  fig_height = fig_width / aspect_ratio
  fig_size = [fig_width, fig_height]

  plt.rcParams.update({
    'backend': 'ps',
    'lines.linewidth': 0.5,
    'font.size': labelsize_pt,
    'axes.titlesize': labelsize_pt,
    'axes.labelsize': labelsize_pt,
    'legend.fontsize': ticksize_pt,
    'xtick.labelsize': ticksize_pt,
    'ytick.labelsize': ticksize_pt,
    'text.usetex': True,
    'font.family': 'serif',
    'font.serif': 'Computer Modern Roman',
    'figure.figsize': fig_size
  })
