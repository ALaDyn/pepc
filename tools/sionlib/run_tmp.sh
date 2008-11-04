#bin/bash
set -x

export SION_DEBUG=SION_pepc.log
export SION_DEBUG_MASK=2047
export SION_DEBUG_RANK1=0
export SION_DEBUG_RANK2=3

WORKDIR=./tmp

gmake clean all 


cd $WORKDIR

~/src/sionlib/examples/pepc/sionpepc -q 10000 -g -n 4 -p 35799040 parts_dump.006500 parts_dump.006500.sion

export SION_DEBUG=SION_dump.log

~/src/sionlib/install/sionlib_jump/bin/siondump -a parts_dump.006500.sion 


~/src/sionlib/examples/pepc/sionpepc -a -v parts_dump.006500.sion parts_dump.006500.ascii
