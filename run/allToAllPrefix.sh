#!/bin/bash
module load mpi/openmpi/3.1

executable="../build/src/tests/distributed_sorter"
numOfStrings=1000000
numOfIterations=4
byteEncoder=5
generator=1
stringLength=500

for dToNRatio in 1.0 0.8 0.6 0.4 0.2 0.0
do
	for byteEncoder in 5 1
	do
		mpirun --bind-to core --map-by core -report-bindings $executable --size $numOfStrings --numberOfIterations $numOfIterations --byteEncoder $byteEncoder --generator $generator --dToNRatio $dToNRatio --stringLength $stringLength --strongScaling
	done
done

#mpirun -np 2 $executable --size $numOfStrings --numberOfIterations $numOfIterations --byteEncoder $byteEncoder --generator $generator --dToNRatio $dToNRatio --stringLength $stringLength
#
#dToNRatio=0.4
#mpirun -np 2 $executable --size $numOfStrings --numberOfIterations $numOfIterations --byteEncoder $byteEncoder --generator $generator --dToNRatio $dToNRatio --stringLength $stringLength
#
#dToNRatio=0.8
#mpirun -np 2 $executable --size $numOfStrings --numberOfIterations $numOfIterations --byteEncoder $byteEncoder --generator $generator --dToNRatio $dToNRatio --stringLength $stringLength
#
#dToNRatio=1.0
#mpirun -np 2 $executable --size $numOfStrings --numberOfIterations $numOfIterations --byteEncoder $byteEncoder --generator $generator --dToNRatio $dToNRatio --stringLength $stringLength
#
#byteEncoder=1
#dToNRatio=0.2
#mpirun -np 2 $executable --size $numOfStrings --numberOfIterations $numOfIterations --byteEncoder $byteEncoder --generator $generator --dToNRatio $dToNRatio --stringLength $stringLength
#
#dToNRatio=0.4
#mpirun -np 2 $executable --size $numOfStrings --numberOfIterations $numOfIterations --byteEncoder $byteEncoder --generator $generator --dToNRatio $dToNRatio --stringLength $stringLength
#
#dToNRatio=0.8
#mpirun -np 2 $executable --size $numOfStrings --numberOfIterations $numOfIterations --byteEncoder $byteEncoder --generator $generator --dToNRatio $dToNRatio --stringLength $stringLength
#
#dToNRatio=1.0
#mpirun -np 2 $executable --size $numOfStrings --numberOfIterations $numOfIterations --byteEncoder $byteEncoder --generator $generator --dToNRatio $dToNRatio --stringLength $stringLength


