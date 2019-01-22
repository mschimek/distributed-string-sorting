#!/bin/bash

touch jobIds.txt
rm jobIds.txt
touch jobIds.txt

sbatch --partition singlenode --nodes=1 --ntasks-per-node=2 --time=00:20:00 $1  >> jobIds.txt
sbatch --partition singlenode --nodes=1 --ntasks-per-node=4 --time=00:20:00 $1  >> jobIds.txt
sbatch --partition singlenode --nodes=1 --ntasks-per-node=8 --time=00:20:00 $1  >> jobIds.txt
sbatch --partition singlenode --nodes=1 --ntasks-per-node=16 --time=00:20:00 $1 >> jobIds.txt
