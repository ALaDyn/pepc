#!/bin/bash
# a script to remove the first (non-balanced) timestep from PEPC's output
#
awk '
 BEGIN { prtln = 1; fl = 1;}
 /.*computing step\s*:\s*0\s*/ { print "FOUND STEP 0, SKIPPING"; prtln = 0; }
 /.*computing step\s*:\s*1\s*/ { print "FOUND STEP 1, ECHOING"; prtln = 1; }
 {if ( prtln == 1 ) { print $0; }}
' $1 > ${1}.filtered
