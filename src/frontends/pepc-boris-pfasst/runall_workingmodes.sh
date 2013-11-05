#!/bin/bash

if [ $# -ne 1 ]
then
  echo "Usage: `basename $0` PARAMFILE"
  exit -1
fi

for i in 1 2 3;
do
  echo "################ $i #################"
  sed -i "s/workingmode[ *]=[ *][1234567890*]/workingmode = $i/" $1
  grep workingmode $1
  ./pepc-boris-pfasst $1
done
