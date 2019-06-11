duration=3
dir="DToNSkewed_$1"
mkdir $dir
cd $dir
../startjobs_multinode_full.sh ../prefix_doubling_skewed.sh
sleep $duration
mv jobIdsMultinodeFull.txt jobIdsPDSkewed_GolombEncoding_01

../startjobs_multinode_full.sh ../merge_sort_skewed.sh
sleep $duration
mv jobIdsMultinodeFull.txt jobIdsMSSkewed_ByteEncoder15_Sampling_23

../startjobs_multinode_full.sh ../kurpicz_skewed.sh
sleep $duration
mv jobIdsMultinodeFull.txt jobIdsKurpiczSkewed

../startjobs_multinode_full.sh ../hQuicksort_skewed.sh
sleep $duration
mv jobIdsMultinodeFull.txt jobIdsHQuicksortSkedwed
