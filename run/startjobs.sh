#!/bin/bash
touch jobIdsMultinode.txt
rm jobIdsMultinode.txt
touch jobIdsMultinode.txt

#sbatch --partition multinode --nodes=2 --ntasks-per-node=1 --time=00:30:00 $1 >> jobIdsMultinode.txt
#sbatch --partition multinode --nodes=4 --ntasks-per-node=1 --time=00:30:00 $1 >> jobIdsMultinode.txt
sbatch --partition multinode --nodes=8 --ntasks-per-node=1 --time=00:30:00 $1 >> jobIdsMultinode.txt
#sbatch --partition multinode --nodes=16 --ntasks-per-node=1 --time=00:30:00 $1 >> jobIdsMultinode.txt

