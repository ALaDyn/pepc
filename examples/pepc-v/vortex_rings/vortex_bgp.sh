#!/bin/sh

if [ $# -ne 4 ]; then
    echo 1>&2 "ERROR: use" $0 "DIR EXEC INPUT SUBMIT"
    exit 127
fi

DIR=$1
EXEC=$2
INPUT=$3
SUBM=$4

if [ -d "$DIR" ]; then
    echo "WARNING:" $DIR "already exists! Are you sure (yes/no)?"
    read answ
    if [ "$answ" != "yes" ]; then
	echo 1>&2 "Aborting..."
	exit 126
    fi
fi

echo 'Preparing' $DIR '...'
if [ ! -d $DIR ]; then
    mkdir $DIR
    mkdir $DIR/log
    mkdir $DIR/part_data
    mkdir $DIR/diag
fi

if [ ! -d $DIR ]; then
    echo 1>&2 "ERROR:" $DIR "was not created, exiting..."
    exit 125
fi

if [ ! -e "$EXEC" ]; then
    echo 1>&2 "ERROR: executable" $EXEC "is missing, exiting..."
    exit 125
fi

if [ ! -e "$INPUT" ]; then
    echo 1>&2 "ERROR: input file" $INPUT "is missing, exiting..."
    exit 124
fi

if [ ! -e "$SUBM" ]; then
    echo 1>&2 "ERROR:" $SUBM "must exist, exiting..."
    exit 123
fi

echo 'Getting executable' $EXEC '...'
cp $EXEC $DIR/pepc
echo 'Getting input file' $INPUT '...'
cp $INPUT $DIR/run.h
echo 'Getting job script' $SUBM '...'
cp $SUBM $DIR/.
echo 'Submitting job in' $DIR '...'
cd $DIR
llsubmit $SUBM
exit 0
