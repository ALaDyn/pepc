export ELG_BLACKLIST=blist
echo "Starting KOJAK analysis ..."
cp tune_old.h run.h
llrun -p4 ../src_tune/pepc 
echo "... done"
