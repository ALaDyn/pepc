#!/bin/bash

FIELD_NPARTS=3
FIELD_PE=1
FIELD_FETCHES=9
FIELD_SHIPS=10
FIELD_INTERACTIONS=11
FIELD_MACEVALS=12
FIELD_RELWORK=13

FIELDS="${FIELD_NPARTS} ${FIELD_FETCHES} ${FIELD_SHIPS} ${FIELD_INTERACTIONS} ${FIELD_MACEVALS} ${FIELD_RELWORK}"

for FIELD in ${FIELDS}
do
  FILEOUT="extr.stats.FIELD${FIELD}.dat"
  rm -f ${FILEOUT}
done


for FILENAME in `ls stats.??????`
do
  echo ${FILENAME}

  for FIELD in ${FIELDS}
  do
    FILEOUT="extr.stats.FIELD${FIELD}.dat"
    tail -n 1024 ${FILENAME} | tr -s ' ' | cut -d ' ' -f ${FIELD} | tr '\n' ' ' >> ${FILEOUT}
    echo "" >> ${FILEOUT}
  done
done

echo "finished."
