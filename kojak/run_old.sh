export ELG_BLACKLIST=blist
echo "Starting KOJAK analysis ..."
cp tune_old.h run.h
llrun -p1 ../src_tune/pepc 
echo "... done"
