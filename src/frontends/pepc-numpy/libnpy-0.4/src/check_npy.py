#!/usr/bin/env python
import numpy as np

def disp_array(arrname):
    arr = np.load(arrname+'.npy')
    print "\nThe %s array %s:" % ( arr.dtype, arrname )
    print arr
  
for arrname in ( 'a',  'b',  'c',  'd',  'e' ):
    disp_array(arrname)
