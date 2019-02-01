#!/bin/bash
# This script gathers time-step-resolved data from the log files that should have been created
# during the JUBE benchmark runs. We look for 'time in step', 'tree walk time', and communication
# timers within the tree walk.

if [ "$#" -ne "2" ]; then
   echo -e "\nmissing prameters"
   echo -e "  usage: $0 <path_to_jube_benchmarks> <jube_benchmark_id>\n"
   exit 1
fi

# define timers we look for and their (file)names
timers[1]="time in step"
timers[2]="tree walk time"
timers[3]="tree comm reqs"
timers[4]="tree comm recv"
short_timers[1]='tis'
short_timers[2]='twt'
short_timers[3]='tcq'
short_timers[4]='tcr'
no_timers=4

# create TMP directory to hold all data (TODO: possibly replace by shared data from JUBE?)
mkdir tmp_steps

# construct proper JUBE id w/ leading zeros
longid=`printf %06d $2`

mpi_case_id=1
# find runs to gather data
for rundir in `find $1/$longid/ -name \*run`
do
   # identify MPI case (we use the job script that should have a single identiable line)
   ### mpi_case=`grep intel $rundir/work/job_jureca.slurm | head -n 1 | awk -F '"' '{print $2}'` # <---- this works for 'multi-MPI-scripts'
   mpi_case=`grep -A 5 'Currently Loaded Modules' $rundir/work/[0-9]*.stderr | awk '{col=1; do {if ($col ~ /.*[mM][Vv]?[Aa]?[pP][iI].*/) {print $col}; col++} while (col<=NF)}' | sed -e 's/\//-/g' | head -n 1 ` # <---- this should work for scripts with 'module list'
   mpi_case="${mpi_case:-"default"}"
   echo "$mpi_case"
   # if we have an output file, look for the data
   if test -f $rundir/work/log; then
      # scan all timers
      for timer in `seq 1 $no_timers`
      do
	 # construct better search pattern for grep (the MPI cases have similar base names...)
	 pattern=`echo ${timers[$timer]} | sed 's/ /./g'`
	 grep "$pattern" $rundir/work/log | awk -F ' : ' '{print $2}' > tmp_steps/${short_timers[$timer]}_${mpi_case}__.${mpi_case_id}.dat
	 # check if file is empty, if so delete
	 test -s tmp_steps/${short_timers[$timer]}_${mpi_case}__.${mpi_case_id}.dat || rm tmp_steps/${short_timers[$timer]}_${mpi_case}__.${mpi_case_id}.dat
      done
      # add case to list of cases we want to plot since we have data
      if ! `echo ${mpi_cases[@]} | grep -q "\<$mpi_case\>"`; then
	 mpi_cases[$mpi_case_id]=$mpi_case
      fi
   fi
   mpi_case_id=$(($mpi_case_id+1))
done

# now collect all timings for one MPI case into a single file
for mpi_case in ${mpi_cases[@]}
do
   # keep timers separate though
   for timer in `seq 1 $no_timers`
   do
      #cat tmp_steps/${short_timers[$timer]}_${mpi_case}__.*.dat > tmp_steps/${short_timers[$timer]}_${mpi_case}__.alldat
      for fil in tmp_steps/${short_timers[$timer]}_${mpi_case}__.*.dat
      do
	tail -n +2 $fil >> tmp_steps/${short_timers[$timer]}_${mpi_case}__.alldat
      done
   done
done

# construct plot file
echo "set term pdf color" > tmp_steps/gnu.plot
echo "set out 'timesteps${longid}.pdf'" >> tmp_steps/gnu.plot

for timer in `seq 1 $no_timers`
do
   mpi_case_id=1

   # start by plotting all data combined
   echo "set title '${timers[$timer]} (rank 0)'" >> tmp_steps/gnu.plot
   echo "set ylabel 'time [s]'" >> tmp_steps/gnu.plot

   # plot violin plots of stats

   # plot time series
   echo "set xlabel 'compiler/MPI combinations'" >> tmp_steps/gnu.plot
   echo -n "set xtics rotate by -15 (" >> tmp_steps/gnu.plot
   id=1
   for mpi_case in ${mpi_cases[@]}
   do
      gnuplot_case=`echo $mpi_case | sed 's/_/ + /g'`
      echo -n "'$gnuplot_case' $id, " >> tmp_steps/gnu.plot
      id=$(($id+1))
   done
   echo -n "'nix' 100) " >> tmp_steps/gnu.plot
   echo " nomirror" >> tmp_steps/gnu.plot
   echo "set nokey" >> tmp_steps/gnu.plot
   echo "set xr [*:*]" >> tmp_steps/gnu.plot
   echo "set yr [0:*]" >> tmp_steps/gnu.plot
   id=1
   for mpi_case in ${mpi_cases[@]}
   do
      file="tmp_steps/${short_timers[$timer]}_${mpi_case}__.alldat"
      len=`wc -l $file | awk '{print $1}'`
      echo "set table \$kdensity${id}" >> tmp_steps/gnu.plot
      echo "plot '$file' using (\$1):(1./${len}) smooth kdensity bandwidth .3 w filledcurve y=0" >> tmp_steps/gnu.plot
      #echo "plot '$file' using (\$1):(1./${len}) smooth kdensity w filledcurve y=0" >> tmp_steps/gnu.plot
      echo "unset table" >> tmp_steps/gnu.plot
      echo "loc${id}=${id}" >> tmp_steps/gnu.plot
      id=$(($id+1))
   done
   echo "plot \\" >> tmp_steps/gnu.plot
   id=1
   for mpi_case in ${mpi_cases[@]}
   do
      file="tmp_steps/${short_timers[$timer]}_${mpi_case}__.alldat"
      gnuplot_case=`echo $mpi_case | sed 's/_/ + /g'`
      #echo "\$kdensity${id} u ((\$2/2.)+loc${id}):1 w filledcurve x=loc${id} lt ${id} t '$mpi_case', \$kdensity${id} u (-(\$2/2.)+loc${id}):1 w filledcurve x=loc${id} lt ${id} notitle, '$file' using (loc${id}):1:(.1) w boxplot fillcolor 'white' lc black lw 3 notitle, \\" >> tmp_steps/gnu.plot
      echo "\$kdensity${id} u ((\$2/2.)+loc${id}):1 w filledcurve x=loc${id} lt ${id} t '$gnuplot_case', \$kdensity${id} u (-(\$2/2.)+loc${id}):1 w filledcurve x=loc${id} lt ${id} notitle, '$file' using (loc${id}):1:(.1) w boxplot lc 'black' lw 2 ps 0.2 notitle, \\" >> tmp_steps/gnu.plot
      id=$(($id+1))
   done
   echo "-1 notitle" >> tmp_steps/gnu.plot

   # plot time series
   echo "set xr [0:*]" >> tmp_steps/gnu.plot
   echo "set yr [0:]" >> tmp_steps/gnu.plot
   echo "set key top" >> tmp_steps/gnu.plot
   echo "unset xtics" >> tmp_steps/gnu.plot
   echo "set xtics mirror" >> tmp_steps/gnu.plot
   echo "set xlabel 'timesteps'" >> tmp_steps/gnu.plot
   echo "plot \\" >> tmp_steps/gnu.plot
   for mpi_case in ${mpi_cases[@]}
   do
      id=1
      gnuplot_case=`echo $mpi_case | sed 's/_/ + /g'`
      for file in tmp_steps/${short_timers[$timer]}_${mpi_case}__.*.dat
      do
	 if [ "$id" -eq "1" ]
	 then
	    echo "'$file' u :1 w lp lt $mpi_case_id t '$gnuplot_case', \\" >> tmp_steps/gnu.plot
	 else
	    echo "'$file' u :1 w lp lt $mpi_case_id notitle, \\" >> tmp_steps/gnu.plot
	 fi
	 id=$(($id+1))
      done
      mpi_case_id=$(($mpi_case_id+1))
   done
   echo "0 notitle" >> tmp_steps/gnu.plot

   # plot each MPI case on single page
   mpi_case_id=1
   for mpi_case in ${mpi_cases[@]}
   do
      echo "plot \\" >> tmp_steps/gnu.plot
      id=1
      gnuplot_case=`echo $mpi_case | sed 's/_/ + /g'`
      for file in tmp_steps/${short_timers[$timer]}_${mpi_case}__.*.dat
      do
	 if [ "$id" -eq "1" ]
	 then
	    echo "'$file' u :1 w lp lt $mpi_case_id t '$gnuplot_case', \\" >> tmp_steps/gnu.plot
	 else
	    echo "'$file' u :1 w lp lt $mpi_case_id notitle, \\" >> tmp_steps/gnu.plot
	 fi
	 id=$(($id+1))
      done
      mpi_case_id=$(($mpi_case_id+1))
      echo "0 notitle" >> tmp_steps/gnu.plot
   done
done

# plot the results
gnuplot tmp_steps/gnu.plot

# clear all data
rm -rf tmp_steps
