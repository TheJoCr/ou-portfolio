#!/bin/sh

module load OpenMPI

#refresh the makefile if necessary
cd `dirname $0`/debug
cmake ..

#Then switch back to the higher directory and build
cd ..
rm ./output/*
cmake --build ./debug
sbatch ./Homework2.sbatch
