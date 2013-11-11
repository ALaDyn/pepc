#!/bin/bash

if [ $# -ne 1 ]
then
  echo "Usage: `basename $0` PARAMFILE"
  exit -1
fi

TEMPFILE="/tmp/runall_workingmodes.$RANDOM"

for i in 1 2 3 4 5;
do
  echo "################ $i #################"
  sed "s/workingmode[ *]=[ *][1234567890*]/workingmode = $i/" $1 > $TEMPFILE
  grep workingmode $TEMPFILE
  ./pepc-boris-pfasst $TEMPFILE
  rm $TEMPFILE
done
