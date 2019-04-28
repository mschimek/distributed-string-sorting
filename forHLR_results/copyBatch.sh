#! /bin/bash

scp -r gw1960@forhlr1.scc.kit.edu:/home/fh1-project-kalb/gw1960/distributed-string-sorting/run_forHLR1/$2 ./$2
for line in $(cat ./$2)
do
  echo $line
 ./copyFromForHLR1.sh $1 $line 
done  
rm ./$2
