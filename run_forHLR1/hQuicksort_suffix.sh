#!/bin/bash
#module load mpi/openmpi/3.1
module load mpi/impi/2018
#export I_MPI_HYDRA_BRANCH_COUNT=-1


path=DummyPath
executable="/home/fh1-project-kalb/gw1960/distributed-string-sorting/build/src/executables/hQuicksort"
size=200000000
stringLength=500
numberOfIterations=11
samplingFactor=2
dToNRatio=0.5
generator=4

mpiexec.hydra -bootstrap slurm $executable --size $size --stringLength $stringLength --numberOfIterations $numberOfIterations --dToNRatio $dToNRatio  --generator $generator --strongScaling --path $path
#mpirun --bind-to core --map-by core  $executable --size $size --stringLength $stringLength --numberOfIterations $numberOfIterations --dToNRatio $dToNRatio  
