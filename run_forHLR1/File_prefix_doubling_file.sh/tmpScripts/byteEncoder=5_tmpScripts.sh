#!/bin/bash
#module load mpi/openmpi/3.1
module load mpi/impi/2018
#export I_MPI_HYDRA_BRANCH_COUNT=-1

path=/work/fh1-project-kalb/gw1960/data/byteEncoder=5
executable="/home/fh1-project-kalb/gw1960/distributed-string-sorting/build/src/executables/prefix_doubling"
numOfStrings=5000000
numOfIterations=3
byteEncoder=5
generator=2
MPIRoutine=2
dToNRatio=0.5
stringLength=1000
for sampler in 2
do
	for golombEncoding in 0 1
	do
		#mpirun --mca coll_tuned_use_dynamic_rules 1 --mca coll_tuned_allgatherv_algorithm 1 --bind-to core --map-by core -report-bindings $executable --size $numOfStrings --numberOfIterations $numOfIterations --byteEncoder $byteEncoder --generator $generator --dToNRatio $dToNRatio --stringLength $stringLength --strongScaling
		#mpirun --bind-to core --map-by core $executable --size $numOfStrings --numberOfIterations $numOfIterations --byteEncoder $byteEncoder --generator $generator --dToNRatio $dToNRatio --stringLength $stringLength --strongScaling --golombEncodingPolicy $golombEncoding --path $path --sampleStringsPolicy $sampler --MPIRoutineAllToAll $MPIRoutine
		mpiexec.hydra -bootstrap slurm $executable --size $numOfStrings --numberOfIterations $numOfIterations --byteEncoder $byteEncoder --generator $generator --dToNRatio $dToNRatio --stringLength $stringLength --strongScaling --golombEncodingPolicy $golombEncoding --path $path --sampleStringsPolicy $sampler --MPIRoutineAllToAll $MPIRoutine
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


