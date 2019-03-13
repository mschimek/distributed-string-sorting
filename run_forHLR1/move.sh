#!/bin/bash
dir="$(date +'%Y_%m_%d_H%H_%M')_$2"
echo $dir
mkdir $dir
counter=1

grep -oE '[^ ]+$' $1 > processedJobIds.txt

for jobid in $(cat processedJobIds.txt)
do
  echo "$jobid"
  jobFile="slurm-${jobid}.out"
  newName="nodes_${counter}"
  counter=$(($counter * 2))
  echo "moving ${jobFile} to ${dir}/${newName}"
  mv $jobFile $dir/$newName
done

rm $1
rm processedJobIds.txt
