#!/bin/bash



touch jobIdsMultinodeFull.txt
rm jobIdsMultinodeFull.txt
touch jobIdsMultinodeFull.txt

sbatch --partition singlenode  --ntasks=1    --time=04:25:00 $1  >> jobIdsMultinodeFull.txt
sbatch --partition singlenode  --ntasks=2    --time=03:50:00 $1  >> jobIdsMultinodeFull.txt
sbatch --partition singlenode  --ntasks=4    --time=02:55:00 $1  >> jobIdsMultinodeFull.txt
sbatch --partition singlenode  --ntasks=8    --time=01:45:00 $1  >> jobIdsMultinodeFull.txt
sbatch --partition singlenode  --ntasks=16   --time=00:55:00 $1  >> jobIdsMultinodeFull.txt
sbatch --partition multinode   --ntasks=32   --time=00:45:00 $1  >> jobIdsMultinodeFull.txt
sbatch --partition multinode   --ntasks=64   --time=00:45:00 $1  >> jobIdsMultinodeFull.txt
sbatch --partition multinode   --ntasks=128  --time=00:45:00 $1  >> jobIdsMultinodeFull.txt
