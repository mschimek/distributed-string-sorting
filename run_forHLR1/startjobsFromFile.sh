#!/bin/bash

prefix="\/work\/fh1-project-kalb\/gw1960\/data\/"
#prefix="\.\.\/upload\/"
mkdir jobIds

for inputFile in $(cat $2)
do
	jobIds="jobIds/$1_${inputFile}.txt"
        touch $jobIds
        rm $jobIds
        touch $jobIds
	completeInputFilePath="$prefix$inputFile"
	echo $completeInputFilePath
	outputFile="./tmpScripts/${inputFile}_tmpScripts.sh"
	echo "outputFile:  $outputFile"

	sed s/DummyPath/"$completeInputFilePath"/ $1 > $outputFile

	#sbatch --partition singlenode  --ntasks=1    --time=08:25:00 $outputFile >> $jobIds
	#sbatch --partition singlenode  --ntasks=2    --time=08:40:00 $outputFile >> $jobIds
	#sbatch --partition singlenode  --ntasks=4    --time=08:25:00 $outputFile >> $jobIds
	#sbatch --partition singlenode  --ntasks=10    --time=08:45:00 $outputFile >> $jobIds
	#sbatch --partition singlenode  --ntasks=20   --time=07:40:00 $outputFile >> $jobIds
	#sbatch --partition multinode   --ntasks=40   --time=06:45:00 $outputFile >> $jobIds
	#sbatch --partition multinode   --ntasks=80   --time=00:10:00 $outputFile >> $jobIds
	#sbatch --partition multinode   --ntasks=160  --time=00:15:00 $outputFile >> $jobIds
	sbatch --partition multinode   --ntasks=320  --time=01:45:00 $outputFile >> $jobIds
	#sbatch --partition multinode   --ntasks=640  --time=01:45:00 $outputFile >> $jobIds
	#sbatch --partition multinode   --ntasks=1280  --time=00:55:00 $outputFile >> $jobIds
	#sbatch --partition multinode   --ntasks=2560  --time=00:55:00 $outputFile >> $jobIds
done
