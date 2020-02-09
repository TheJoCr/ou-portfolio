#!/bin/sh
set -e

#refresh the makefile if necessary
cd `dirname $0`/debug
cmake ..

#Then switch back to the higher directory and build
cd ..
rm -f ./output/*
cmake --build ./debug

sbatch Homework5.sbatch