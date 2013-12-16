#!/bin/bash

vtkdir="vtk"

if [ ! -n "$1" ]; then

  echo " "
  echo "Usage: call"
  echo " "
  echo "       pepc_particle2vtk.sh"
  echo " "
  echo "from the directory containing the particle_xxxxx_yyyy.dat files."
  echo " "
  echo "The script will read them, collect the data into a single"
  echo "vtk-file for each timestep and put these into the subdir $vtkdir"
  echo " "
  echo "These files can then be loaded with ParaView for example."
  echo " "



  rm -rf ./$vtkdir
  mkdir ./$vtkdir
  cd $vtkdir

  for timestep in `ls ../particle_*_000000.dat | awk -F _ '{print $2}'` ; do
    ../$0 $timestep
  done
else

  timestep=$1
  echo "==== Processing timestep $timestep ===="
  # merge all n-1 input files
  # check for *.tmp files
  fileout="particles_$timestep.vtk"

  for cpu in `ls ../particle_${timestep}_*.dat | awk -F _ '{print $3}' | awk -F . '{print $1}'` ; do

    filein="../particle_${timestep}_${cpu}.dat"

    # echo "file to merge: ${filein}"

     # merge all files together and btw replace NaN by 0
    awk -v c=$cpu '{ gsub("NaN", "0.99999E30", $0); print c, $0}' < $filein | grep -v "#" >> $fileout.tmp
  done

  numparticles=`wc -l < $fileout.tmp`
  let seqlimit=$numparticles-1
  let numvertices=$numparticles*2

  echo "found $numparticles particles in total"

  echo "# vtk DataFile Version 3.0" > $fileout
  echo "vtk output" >> $fileout
  echo "ASCII" >> $fileout

  # Read particle coordinates as "grid-points"
  echo "DATASET POLYDATA" >> $fileout
  echo "POINTS $numparticles float" >> $fileout
  awk '{print $3 " " $4 " " $5}' < $fileout.tmp >> $fileout

  # list of vertices, i.e. lines with "1 xxx" in them, xxx running from 0 to number of points-1
  echo "" >> $fileout
  echo "VERTICES $numparticles $numvertices" >> $fileout
  echo -e " 1 \c" >> $fileout
  echo -e `seq -s "\n 1 " 0 $seqlimit` >> $fileout

  echo "" >> $fileout
  echo "POINT_DATA $numparticles" >> $fileout


  # particle charge
  echo "" >> $fileout
  echo "SCALARS q float" >> $fileout
  echo "LOOKUP_TABLE default" >> $fileout
  awk '{print $9}' < $fileout.tmp >> $fileout

  # particle mass
  echo "" >> $fileout
  echo "SCALARS m float" >> $fileout
  echo "LOOKUP_TABLE default" >> $fileout
  awk '{print $10}' < $fileout.tmp >> $fileout

  # Read particle velocities
  echo "" >> $fileout
  echo "SCALARS u_velocity float 3" >> $fileout
  echo "LOOKUP_TABLE default" >> $fileout
  awk '{print $6 " " $7 " " $8}' < $fileout.tmp >> $fileout



 # rm -f $fileout.tmp

fi
