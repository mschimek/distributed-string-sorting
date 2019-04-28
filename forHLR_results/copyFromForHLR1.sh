#! /bin/bash

scp -r gw1960@forhlr1.scc.kit.edu:/home/fh1-project-kalb/gw1960/distributed-string-sorting/run_forHLR1/${2} ./forHLR1/
./processDataToInputFile.sh forHLR1/${2}
Rscript ${1} forHLR1/${2}
