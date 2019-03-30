#!/bin/bash
module load mpi/openmpi/3.1

path=DummyPath
executable="../build/src/executables/prefix_doubling"
numOfStrings=500000
numOfIterations=10
byteEncoder=5
generator=2
MPIRoutine=2
dToNRatio=0.5
stringLength=1000

for sampler in 1 0
do
	for golombEncoding in 0 1 2
	do
		#mpirun --mca coll_tuned_use_dynamic_rules 1 --mca coll_tuned_allgatherv_algorithm 1 --bind-to core --map-by core -report-bindings $executable --size $numOfStrings --numberOfIterations $numOfIterations --byteEncoder $byteEncoder --generator $generator --dToNRatio $dToNRatio --stringLength $stringLength --strongScaling
		mpirun --bind-to core --map-by core $executable --size $numOfStrings --numberOfIterations $numOfIterations --byteEncoder $byteEncoder --generator $generator --dToNRatio $dToNRatio --stringLength $stringLength --strongScaling --golombEncodingPolicy $golombEncoding --path $path --sampleStringsPolicy $sampler --MPIRoutineAllToAll $MPIRoutine
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

