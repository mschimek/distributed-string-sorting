#!/bin/bash

prefix="\/work\/fh1-project-kalb\/gw1960\/data\/"
#prefix="\.\.\/upload\/"
dir="./jobIds"

mkdir $dir

for inputFile in $(cat $2)
do
	jobIds="./jobIds/$1_${inputFile}.txt"
	echo "jobIds: $jobIds"
        touch $jobIds
        rm $jobIds
        touch $jobIds
	completeInputFilePath="$prefix$inputFile"
	echo $completeInputFilePath
	outputFile="./tmpScripts/${inputFile}_tmpScripts.sh"
	echo "outputFile:  $outputFile"
	
	executable="../$1"

	sed s/DummyPath/"$completeInputFilePath"/ $executable > $outputFile

	sbatch --partition develop  --ntasks=2   --time=00:05:00 $outputFile >> $jobIds
	#sbatch --partition singlenode  --ntasks=2    --time=08:40:00 $outputFile >> $jobIds
	#sbatch --partition singlenode  --ntasks=4    --time=08:25:00 $outputFile >> $jobIds
	#sbatch --partition singlenode  --ntasks=10    --time=08:45:00 $outputFile >> $jobIds
	#sbatch --partition develop  --ntasks=20   --time=00:04:00 $outputFile >> $jobIds
	#sbatch --partition multinode   --ntasks=40   --time=0:05:00 $outputFile >> $jobIds
	#sbatch --partition multinode   --ntasks=80   --time=00:10:00 $outputFile >> $jobIds
	#sbatch --partition multinode   --ntasks=160  --time=01:45:00 $outputFile >> $jobIds
	#sbatch --partition multinode   --ntasks=320  --time=01:45:00 $outputFile >> $jobIds
	#sbatch --partition multinode   --ntasks=480  --time=01:45:00 $outputFile >> $jobIds
	#sbatch --partition multinode   --ntasks=640  --time=01:45:00 $outputFile >> $jobIds
	#sbatch --partition multinode   --ntasks=960  --time=01:45:00 $outputFile >> $jobIds
	#sbatch --partition multinode   --ntasks=1280  --time=01:55:00 $outputFile >> $jobIds
	#sbatch --partition multinode   --ntasks=2560  --time=00:55:00 $outputFile >> $jobIds
done
