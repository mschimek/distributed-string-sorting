#!/bin/bash
dir="$(date +'%Y_%m_%d_H%H_%M')_$1"
echo $dir
mkdir $dir
counter=1

grep -oE '[^ ]+$' jobIds.txt > processedJobIds.txt

for jobid in $(cat processedJobIds.txt)
do
  echo "$jobid"
  counter=$(($counter * 2))
  jobFile="slurm-${jobid}.out"
  newName="nodes_${counter}"
  echo "moving ${jobFile} to ${dir}/${newName}"
  mv $jobFile $dir/$newName
done

rm jobIds.txt
rm processedJobIds.txt
