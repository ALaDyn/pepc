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

firststats=`ls stats.* | head -n 1`

NUMPE=`grep "# procs, walk_threads, max_particles_per_thread:" stats.000000 | tr -s ' ' | cut -d ' ' -f 6`

echo "Found data for ${NUMPE} processors"

for FILENAME in `ls stats.??????`
do
  echo ${FILENAME}

  for FIELD in ${FIELDS}
  do
    FILEOUT="extr.stats.FIELD${FIELD}.dat"
    tail -n ${NUMPE} ${FILENAME} | tr -s ' ' | cut -d ' ' -f ${FIELD} | tr '\n' ' ' >> ${FILEOUT}
    echo "" >> ${FILEOUT}
  done
done

echo "finished."
