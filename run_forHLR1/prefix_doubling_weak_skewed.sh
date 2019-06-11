#!/bin/bash
#module load mpi/openmpi/3.1
module load mpi/impi/2018
#export I_MPI_HYDRA_BRANCH_COUNT=-1

executable="/home/fh1-project-kalb/gw1960/distributed-string-sorting/build/src/executables/prefix_doubling"
numOfStrings=500000
numOfIterations=6
sampler=2
byteEncoder=5
MPIRoutine=2
generator=3
stringLength=500

for dToNRatio in 0.0 0.25 0.5 0.75 1.0
do
	for golombEncoding in 0 1
	do
		#mpirun --mca coll_tuned_use_dynamic_rules 1 --mca coll_tuned_allgatherv_algorithm 1 --bind-to core --map-by core -report-bindings $executable --size $numOfStrings --numberOfIterations $numOfIterations --byteEncoder $byteEncoder --generator $generator --dToNRatio $dToNRatio --stringLength $stringLength --strongScaling
		#mpirun --bind-to core --map-by core $executable --size $numOfStrings --numberOfIterations $numOfIterations --byteEncoder $byteEncoder --generator $generator --dToNRatio $dToNRatio --stringLength $stringLength --strongScaling --golombEncodingPolicy $golombEncoding --sampleStringsPolicy $sampler --MPIRoutineAllToAll $MPIRoutine
		mpiexec.hydra -bootstrap slurm $executable --size $numOfStrings --numberOfIterations $numOfIterations --byteEncoder $byteEncoder --generator $generator --dToNRatio $dToNRatio --stringLength $stringLength --golombEncodingPolicy $golombEncoding --sampleStringsPolicy $sampler --MPIRoutineAllToAll $MPIRoutine --compressLcps
	done
done
