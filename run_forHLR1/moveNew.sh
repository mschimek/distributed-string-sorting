#!/bin/bash
ls >> outputTMP4848
grep  -rE '^jobIds' outputTMP4848 > ListJobIds
rm outputTMP4848
touch writtenResults
for jobIds in $(cat ListJobIds)
do
  counter=1
  dir="${jobIds#"jobIds"}"
  echo $dir >> writtenResults
  mkdir $dir
  grep -oE '[^ ]+$' $jobIds > processedJobIds.txt
  for jobid in $(cat processedJobIds.txt)
  do
    jobFile="slurm-${jobid}.out"
    newName="nodes_${counter}"
    counter=$(($counter * 2))
    mv $jobFile $dir/$newName
  done
  rm processedJobIds.txt
done
rm ListJobIds
