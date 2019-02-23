#!/bin/bash
module load mpi/impi/2019

executable="../build/src/tests/mpitest"
sizeAllGather=10000
sizeAllToAll=500000000
numberOfIterations=20

for ratio in {10..0..2}
do
	#mpirun --mca coll_tuned_use_dynamic_rules 1 --mca coll_tuned_allgatherv_algorithm 1 --bind-to core --map-by core -report-bindings $executable --size $numOfStrings --numberOfIterations $numOfIterations --byteEncoder $byteEncoder --generator $generator --dToNRatio $dToNRatio --stringLength $stringLength --strongScaling
	actualSizeAllToAll=$(($ratio * $sizeAllToAll / 10))
	mpirun --bind-to core --map-by core $executable --sizeAllGather $sizeAllGather --sizeAllToAll $actualSizeAllToAll --numberOfIterations $numberOfIterations
done

