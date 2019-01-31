#!/bin/bash

touch jobIdsSingleNode.txt
rm jobIdsSingleNode.txt
touch jobIdsSingleNode.txt

sbatch --partition singlenode  --ntasks=2    --time=00:20:00 $1  >> jobIdsSingleNode.txt
sbatch --partition singlenode  --ntasks=4    --time=00:20:00 $1  >> jobIdsSingleNode.txt
sbatch --partition singlenode  --ntasks=8    --time=00:20:00 $1  >> jobIdsSingleNode.txt
sbatch --partition singlenode  --ntasks=16   --time=00:20:00 $1  >> jobIdsSingleNode.txt
