#!/bin/bash
module load mpi/mvapich2/2.3

executable="../build/src/tests/mpitest"

for dToNRatio in 1.0 0.8 0.6 0.4 0.2 0.0
do
	for byteEncoder in 5 1
	do
		mpiexec $executable 
		echo "done"
	done
done

