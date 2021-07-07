#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 30 09:30:12 2018

@author: hawkes
"""

import os
import sys
import fileinput

variable_file = 'update_variable.txt'
total_iteration = int(sys.argv[1])

#===========================Read-in and record updated values================
variable_f = fileinput.input(variable_file,False)

var_values = []
for line in variable_f:
    var_values.append(line.rstrip())

fileinput.close()

#=========================Read-in templates and open new file================
input_file_name = './params_template'
output_file_name = 'params'

ofile = open(output_file_name, 'w')
f = fileinput.input(input_file_name, False)

#===================Write previously read variables into new file============
remaining_step = total_iteration - int(var_values[1])
#if (remaining_step > 1000000):
#    remaining_step = 1000000
#print(remaining_step)
checkW = ('REPLACE_RESUME', 'REPLACE_TIMEIN', 'REPLACE_TNP', 'REPLACE_VP1', 'REPLACE_VP2', 'REPLACE_NSTEP')
replaceW = (var_values[0], var_values[1], var_values[2], var_values[3], var_values[4], str(remaining_step))

for line in f:
    for check, rep in zip(checkW, replaceW):
        line = line.replace(check, rep)
    ofile.write(line)

#sys.stderr.write(var_values[2])
ofile.close()
fileinput.close()

#==================Perform the same with job_template========================
#input_file_name = 'job_template'
#output_file_name = 'job_script'
#
#ofile = open(output_file_name, 'w')
#f = fileinput.input(input_file_name, False)
#
#nodes = 8
#
#checkW = ('N')
#replaceW = (str(nodes))
#
#for line in f:
#    for check, rep in zip(checkW, replaceW):
#        line = line.replace(check, rep)
#    ofile.write(line)
#
#ofile.close()
#fileinput.close()
