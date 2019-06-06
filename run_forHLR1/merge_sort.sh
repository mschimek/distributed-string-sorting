#!/bin/bash
module load mpi/openmpi/3.1
#module load mpi/impi/2018
#export I_MPI_HYDRA_BRANCH_COUNT=-1

executable="../build/src/executables/distributed_sorter"
numOfStrings=200000000
sampler=2
numOfIterations=5
byteEncoder=1
generator=1
stringLength=500
MPIRoutine=2
dToNRatio=0.5

for dToNRatio in 0.0 0.25 0.5 0.75 1.0
do
	for byteEncoder in 1
	do
		#mpirun --mca coll_tuned_use_dynamic_rules 1 --mca coll_tuned_allgatherv_algorithm 1 --bind-to core --map-by core -report-bindings $executable --size $numOfStrings --numberOfIterations $numOfIterations --byteEncoder $byteEncoder --generator $generator --dToNRatio $dToNRatio --stringLength $stringLength --strongScaling
                #mpiexec.hydra -bootstrap slurm  $executable --size $numOfStrings --numberOfIterations $numOfIterations --byteEncoder $byteEncoder --generator $generator --dToNRatio $dToNRatio --stringLength $stringLength --strongScaling  --sampleStringsPolicy $sampler --MPIRoutineAllToAll $MPIRoutine --compressLcps
		mpirun --bind-to core --map-by core  $executable --size $numOfStrings --numberOfIterations $numOfIterations --byteEncoder $byteEncoder --generator $generator --dToNRatio $dToNRatio --stringLength $stringLength --strongScaling  --sampleStringsPolicy $sampler --MPIRoutineAllToAll $MPIRoutine --compressLcps
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


