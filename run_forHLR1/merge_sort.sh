#!/bin/bash
module load mpi/openmpi/3.1

executable="../build/src/executables/distributed_sorter"
numOfStrings=5000000
sampler=0
numOfIterations=10
byteEncoder=0
generator=1
stringLength=1000
MPIRoutine=2

for dToNRatio in 0.0 0.25 0.5 0.75 1.0
do
	for byteEncoder in 1 5
	do
		#mpirun --mca coll_tuned_use_dynamic_rules 1 --mca coll_tuned_allgatherv_algorithm 1 --bind-to core --map-by core -report-bindings $executable --size $numOfStrings --numberOfIterations $numOfIterations --byteEncoder $byteEncoder --generator $generator --dToNRatio $dToNRatio --stringLength $stringLength --strongScaling
		mpirun --bind-to core --map-by core  $executable --size $numOfStrings --numberOfIterations $numOfIterations --byteEncoder $byteEncoder --generator $generator --dToNRatio $dToNRatio --stringLength $stringLength --strongScaling  --sampleStringsPolicy $sampler --MPIRoutineAllToAll $MPIRoutine
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


