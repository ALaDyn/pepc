#!/bin/sh
export ELG_BLACKLIST=blist
echo "Starting KOJAK analysis ..."
echo "Using inputs from kojak.h"
cp kojak.h run.h
llrun -p4 ../pepc-b/pepck 
echo "... done"
