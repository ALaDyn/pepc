#!/bin/bash
#
for d in `find $1 -name job.ompi.slurm.skel | xargs dirname`;
do
   grep -H 'bound to' $d/[1-9]*.stderr | sed -e 's/\[/:/' -e 's/W /W:/' -e 's/ bo/:bo/' -e 's/\.juwels\.fzj\.de//' | sort -k3n,3n
done | awk -F ":" '{print $1, $3, $5,"\t"$7}'
