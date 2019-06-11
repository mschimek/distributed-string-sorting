#!/bin/bash
#module load mpi/openmpi/3.1
module load mpi/impi/2018
#export I_MPI_HYDRA_BRANCH_COUNT=-1


executable="/home/fh1-project-kalb/gw1960/distributed-string-sorting/build/src/executables/hQuicksort"
size=500000
stringLength=500
numberOfIterations=6
samplingFactor=2
dToNRatio=0.5
generator=3

for dToNRatio in 0.0 0.25 0.5 0.75 1.0 
do
		mpiexec.hydra -bootstrap slurm $executable --size $size --stringLength $stringLength --numberOfIterations $numberOfIterations --dToNRatio $dToNRatio  --generator $generator
		#mpirun --bind-to core --map-by core  $executable --size $size --stringLength $stringLength --numberOfIterations $numberOfIterations --dToNRatio $dToNRatio  
done
