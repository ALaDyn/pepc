#!/bin/sh
# for help:
#
#mpirun -h


#
# to run with dynamically allocated partition:
#
TASKS=$1  #  change to desired number of MPI tasks
MODE=$2
TOP=MESH  # Topology - MESH or TORUS
PEPC=../pepc-b/pepcb
cp bench.h run.h
echo $PEPC
mpirun -np $TASKS -connect $TOP -mode $MODE -exe `/bin/pwd`/$PEPC -cwd `/bin/pwd`


#
# to run with a pre-allocated partition:
#
PARTITION=R000_N0  #  change to your partition name
#mpirun -partition $PARTITION -exe `/bin/pwd`/helloworld.rts -cwd `/bin/pwd`

