duration=1
dir="DToNSkewed_$1"
mkdir $dir
cd $dir
#../startjobs_multinode_full.sh ../prefix_doubling_skewed.sh
#sleep $duration
#mv jobIdsMultinodeFull.txt jobIdspd

../startjobs_multinode_full.sh ../merge_sort_skewed.sh
sleep $duration
mv jobIdsMultinodeFull.txt jobIdsms

#../startjobs_multinode_full.sh ../kurpicz_skewed.sh
#sleep $duration
#mv jobIdsMultinodeFull.txt jobIdskurpicz
#
#../startjobs_multinode_full.sh ../hQuicksort_skewed.sh
#sleep $duration
#mv jobIdsMultinodeFull.txt jobIdshquick
