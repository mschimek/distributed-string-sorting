#!/bin/bash
dir="$2"
echo $dir
mkdir $dir
counter=1

grep -oE '[^ ]+$' $1 > processedJobIds.txt

for jobid in $(cat processedJobIds.txt)
do
  jobFile="slurm-${jobid}.out"
  newName="nodes_${counter}"
  counter=$(($counter * 2))
  mv $jobFile $dir/$newName
done

rm $1
rm processedJobIds.txt
