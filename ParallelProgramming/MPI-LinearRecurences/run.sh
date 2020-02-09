#!/bin/sh

module load OpenMPI

#refresh the makefile if necessary
cd `dirname $0`/debug
cmake ..

#Then switch back to the higher directory and build
cd ..
#rm ./output/*
cmake --build ./debug

NTASKS=(1 2 4 8 16 32)
NPERNODE=(1 2 4 8 16 16)
NNODES=(1 1 1 1 1 2)

for (( i=1; i<8; i++ ));
do
  echo "Starting Job" $i " :  ntasks=" ${NTASKS[$i-1]} " per-node=" ${NPERNODE[$i-1]} " nodes=" ${NNODES[$i-1]}
  sbatch --ntasks ${NTASKS[$i-1]} --ntasks-per-node ${NPERNODE[$i-1]} --nodes ${NNODES[$i-1]} ./Homework3.sbatch
done