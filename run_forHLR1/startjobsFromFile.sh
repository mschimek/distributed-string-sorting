#!/bin/bash

prefix="\.\.\/upload\/"

for inputFile in $(cat $2)
do
	jobIds="jobIds_${inputFile}.txt"
        touch $jobIds
        rm $jobIds
        touch $jobIds
	completeInputFilePath="$prefix$inputFile"
	echo $completeInputFilePath
	outputFile="./tmpScripts/${inputFile}_tmpScripts.sh"
	echo "outputFile:  $outputFile"

	sed s/DummyPath/"$completeInputFilePath"/ $1 > $outputFile

	sbatch --partition singlenode  --ntasks=1    --time=02:25:00 $outputFile >> $jobIds
	sbatch --partition singlenode  --ntasks=2    --time=00:20:00 $outputFile >> $jobIds
	sbatch --partition singlenode  --ntasks=4    --time=00:25:00 $outputFile >> $jobIds
	sbatch --partition singlenode  --ntasks=8    --time=00:45:00 $outputFile >> $jobIds
	sbatch --partition singlenode  --ntasks=16   --time=00:40:00 $outputFile >> $jobIds
	sbatch --partition multinode   --ntasks=32   --time=00:45:00 $outputFile >> $jobIds
	sbatch --partition multinode   --ntasks=64   --time=00:45:00 $outputFile >> $jobIds
	sbatch --partition multinode   --ntasks=128  --time=00:45:00 $outputFile >> $jobIds
done
