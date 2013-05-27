#!/bin/bash

if [ $# -ne 1 ]
then
  echo "Usage: `basename $0` RUN|TEST"
  exit 5
fi

if [ "$1" == "RUN" ]
then
    FILEIN='runsims.out'
    FILEOUT='runsims.run'

    echo "= Reading from file ${FILEIN}, writing to ${FILEOUT}"

    rm -f ${FILEOUT}

    IFS=$'\n'
    for LINE in `(cat $FILEIN)`
    do
      [[ "$LINE" =~ ^#.*$ ]] && echo -e "${LINE}" >> ${FILEOUT} && continue

      shortline=`echo "$LINE" | tr -s ' '`
      TEEV=`echo "$shortline" | cut -d ' ' -f 2`
      NENM3=`echo "$shortline" | cut -d ' ' -f 3`
      VOVTH=`echo "$shortline" | cut -d ' ' -f 4`
      WWPL=`echo "$shortline" | cut -d ' ' -f 5`
      NE=`echo "$shortline" | cut -d ' ' -f 6`
      TIMES=`echo "$shortline" | cut -d ' ' -f 7`

      echo "== ${TEEV} ${NENM3} ${VOVTH} ${WWPL} ${NE} ${TIMES}"

      DIRNAME="Te${TEEV}_ne${NENM3}_V0VTHERM${VOVTH}_WWPL${WWPL}_NE${NE}"
      echo "== dirname=${DIRNAME}"

      mkdir -p ${DIRNAME}
      cd ${DIRNAME}

      cp ../runsims.h ./run.h

      sed -i "s/#TE#/${TE}/g" run.h
      sed -i "s/#NE#/${NE}/g" run.h
      sed -i "s/#NENM3#/${NENM3}/g" run.h
      sed -i "s/#VOVTH#/${VOVTH}/g" run.h
      sed -i "s/#WWPL#/${WWPL}/g" run.h

      cp ../pepc-collfreq ./pepc-collfreq
      
      let "WCHOURS=5*$TIMES/11000+1"

      cp ../runsims.juqueen ./runsims.juqueen
      sed -i "s/#JOBNAME#/${DIRNAME}/g" runsims.juqueen
      sed -i "s/#WCHOURS#/${WCHOURS}/g" runsims.juqueen
      
      llsubmit runsims.juqueen

      echo "${TEEV} ${NENM3} ${VOVTH} ${WWPL} ${NE} ${TIMES} ${WCHOURS}" >> ../$FILEOUT

      cd ..  
    done
fi

if [ "$1" == "TEST" ]
then
    FILEIN='runsims.dat'
    FILEOUT='runsims.out'

    echo "= Reading from file ${FILEIN}, writing to ${FILEOUT}"

    rm -f ${FILEOUT}

    IFS=$'\n'
    for LINE in `(cat $FILEIN)`
    do
      [[ "$LINE" =~ ^#.*$ ]] && echo -e "${LINE}" >> ${FILEOUT} && continue

      shortline=`echo "$LINE" | tr -s ' '`
      TEEV=`echo "$shortline" | cut -d ' ' -f 2`
      NENM3=`echo "$shortline" | cut -d ' ' -f 3`
      VOVTH=`echo "$shortline" | cut -d ' ' -f 4`
      WWPL=`echo "$shortline" | cut -d ' ' -f 5`
      NE=`echo "$shortline" | cut -d ' ' -f 6`
      TIMES=`echo "$shortline" | cut -d ' ' -f 7`

      echo "== ${TEEV} ${NENM3} ${VOVTH} ${WWPL} ${NE} ${TIMES}"

      DIRNAME="Te${TEEV}_ne${NENM3}_V0VTHERM${VOVTH}_WWPL${WWPL}_NE${NE}"
      echo "== dirname=${DIRNAME}"

      mkdir -p ${DIRNAME}
      cd ${DIRNAME}

      cp ../runsims.h ./run.h

      sed -i "s/#TE#/${TE}/g" run.h
      sed -i "s/#NE#/${NE}/g" run.h
      sed -i "s/#NENM3#/${NENM3}/g" run.h
      sed -i "s/#VOVTH#/${VOVTH}/g" run.h
      sed -i "s/#WWPL#/${WWPL}/g" run.h

      ln -sf ../pepc-collfreq ./pepc-collfreq

      ./pepc-collfreq ./run.h > run.out & PID=$!

      echo "== PID=${PID}"

      sleep 3s

      kill $PID 2>&1 > /dev/null

      NT=`grep " nt |" run.out | tr -d ' ' | cut -d "|" -f 3`

      echo == NT=${NT}

      echo "${LINE}        ${NT}" >> ../$FILEOUT

      cd ..  
    done

    killall pepc-collfreq 2>&1 > /dev/null
fi

