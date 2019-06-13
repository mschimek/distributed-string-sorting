#! /bin/bash

scp -r gw1960@forhlr1.scc.kit.edu:/home/fh1-project-kalb/gw1960/distributed-string-sorting/run_forHLR1/$2/writtenResults writtenResults
mkdir evaluation/$2
for line in $(cat ./writtenResults)
do
  completePath="$2/$line"
  echo "complete path: $completePath"
  scp -r gw1960@forhlr1.scc.kit.edu:/home/fh1-project-kalb/gw1960/distributed-string-sorting/run_forHLR1/$completePath ./evaluation/$completePath/
  ./processDataToInputFile.sh evaluation/$completePath
done  
rm writtenResults
