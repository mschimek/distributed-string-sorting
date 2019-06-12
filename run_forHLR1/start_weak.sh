#!/bin/bash



touch jobIdsMultinodeFull.txt
rm jobIdsMultinodeFull.txt
touch jobIdsMultinodeFull.txt

#sbatch --partition singlenode  --ntasks=1    --time=01:25:00 $1  >> jobIdsMultinodeFull.txt
#sbatch --partition singlenode  --ntasks=2    --time=03:50:00 $1  >> jobIdsMultinodeFull.txt
#sbatch --partition singlenode  --ntasks=4    --time=03:55:00 $1  >> jobIdsMultinodeFull.txt
#sbatch --partition singlenode  --ntasks=10    --time=01:45:00 $1  >> jobIdsMultinodeFull.txt
#sbatch --partition develop  --ntasks=20   --time=00:35:00 $1  >> jobIdsMultinodeFull.txt
#sbatch --partition multinode   --ntasks=40   --time=00:05:00 $1  >> jobIdsMultinodeFull.txt
#sbatch --partition multinode   --ntasks=80   --time=00:05:00 $1  >> jobIdsMultinodeFull.txt
#sbatch --partition multinode   --ntasks=160  --time=00:05:00 $1  >> jobIdsMultinodeFull.txt
#sbatch --partition multinode   --ntasks=320  --time=00:05:00 $1  >> jobIdsMultinodeFull.txt
#sbatch --partition multinode   --ntasks=640  --time=00:05:00 $1  >> jobIdsMultinodeFull.txt
#sbatch --partition multinode   --ntasks=1280  --time=00:40:00 $1  >> jobIdsMultinodeFull.txt

#batch --partition singlenode  --ntasks=1    --time=00:05:00 $1  >> jobIdsMultinodeFull.txt
#sbatch --partition develop  --ntasks=2    --time=00:05:00 $1  >> jobIdsMultinodeFull.txt
#sbatch --partition singlenode  --ntasks=4    --time=00:05:00 $1  >> jobIdsMultinodeFull.txt
#sbatch --partition singlenode  --ntasks=8    --time=00:05:00 $1  >> jobIdsMultinodeFull.txt
sbatch --partition singlenode  --ntasks=20   --time=00:45:00 $1  >> jobIdsMultinodeFull.txt
#sbatch --partition multinode   --ntasks=40   --time=00:45:00 $1  >> jobIdsMultinodeFull.txt
#sbatch --partition multinode   --ntasks=80   --time=00:45:00 $1  >> jobIdsMultinodeFull.txt
#sbatch --partition multinode   --ntasks=160  --time=00:45:00 $1  >> jobIdsMultinodeFull.txt
#sbatch --partition multinode   --ntasks=320  --time=02:15:00 $1  >> jobIdsMultinodeFull.txt
#sbatch --partition multinode   --ntasks=640  --time=02:15:00 $1  >> jobIdsMultinodeFull.txt
#sbatch --partition multinode   --ntasks=1280  --time=02:35:00 $1  >> jobIdsMultinodeFull.txt
#sbatch --partition multinode   --ntasks=2560  --time=01:35:00 $1  >> jobIdsMultinodeFull.txt
#


#sbatch --partition multinode   --ntasks=320  --time=00:35:00 $1  >> jobIdsMultinodeFull.txt
#sbatch --partition multinode   --ntasks=500  --time=00:15:00 $1  >> jobIdsMultinodeFull.txt
#sbatch --partition multinode   --ntasks=1000  --time=00:15:00 $1  >> jobIdsMultinodeFull.txt
