duration=1
dir="DToN_$1"
mkdir $dir
cd $dir
../startjobs_multinode_full.sh ../prefix_doubling.sh
sleep $duration
mv jobIdsMultinodeFull.txt jobIdspd

../startjobs_multinode_full.sh ../merge_sort.sh
sleep $duration
mv jobIdsMultinodeFull.txt jobIdsms

../startjobs_multinode_full.sh ../kurpicz.sh
sleep $duration
mv jobIdsMultinodeFull.txt jobIdskurpicz

../startjobs_multinode_full.sh ../hQuicksort.sh
sleep $duration
mv jobIdsMultinodeFull.txt jobIdshquick
