#!/bin/bash



prefix="../upload/"

for file in $(cat $2)
do
	echo file
	filename="jobIds_$file.txt"
        touch $filename 
        rm $filename
        touch $filename
	completeFileName="$prefix$filename"
	echo $completeFileName

	sed s/$completeFileName/DummyPath > tmpScript

	echo $(cat tmpScript)
	sbatch --partition singlenode  --ntasks=1    --time=02:25:00 tmpScript  >> $filename
	#sbatch --partition singlenode  --ntasks=2    --time=00:20:00 $1  >> $filename
	#sbatch --partition singlenode  --ntasks=4    --time=00:25:00 $1  >> $filename
	#sbatch --partition singlenode  --ntasks=8    --time=00:45:00 $1  >> $filename
	#sbatch --partition singlenode  --ntasks=16   --time=00:40:00 $1  >> $filename
	#sbatch --partition multinode   --ntasks=32   --time=00:45:00 $1  >> $filename
	#sbatch --partition multinode   --ntasks=64   --time=00:45:00 $1  >> $filename
	#sbatch --partition multinode   --ntasks=128  --time=00:45:00 $1  >> $filename
	rm tmpScript
done
