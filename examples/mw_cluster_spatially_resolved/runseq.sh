#!/bin/bash

if [ $# -ne 2 ]
then
  echo "Usage: `basename $0` START# END#"
  exit 1
fi

for i in `seq $1 $2`; do
  ival=`printf "%5.5d" ${i}`
  ./runsingle.sh ${ival}
done


