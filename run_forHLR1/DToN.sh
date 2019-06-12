duration=3
dir="DToN_$1"
mkdir $dir
cd $dir
../startjobs_multinode_full.sh ../prefix_doubling.sh
sleep $duration
mv jobIdsMultinodeFull.txt jobIdsPD_GolombEncoding_01_200000000_500

../startjobs_multinode_full.sh ../merge_sort.sh
sleep $duration
mv jobIdsMultinodeFull.txt jobIdsMS_ByteEncoder15_Sampling_23_200000000_500

../startjobs_multinode_full.sh ../kurpicz.sh
sleep $duration
mv jobIdsMultinodeFull.txt jobIdsKurpicz_200000000_500

../startjobs_multinode_full.sh ../hQuicksort.sh
sleep $duration
mv jobIdsMultinodeFull.txt jobIdsHQuicksort_200000000_500
