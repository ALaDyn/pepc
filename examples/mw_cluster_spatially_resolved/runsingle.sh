#!/bin/bash

echo $1

cp -r rundir-00000 rundir-$1
cd rundir-$1

touch run.h
echo "&pepcmw"     >> run.h
echo "rngseed=$1"  >> run.h
cat run.h.template >> run.h

msub cluster.juropa
cd ..
