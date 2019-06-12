duration=1
dir="DToNWeak_$1"
mkdir $dir
cd $dir
../start_weak.sh ../prefix_doubling_weak.sh
sleep $duration
mv jobIdsMultinodeFull.txt jobIdsPDWeak_GolombEncoding_01

../start_weak.sh ../merge_sort_weak.sh
sleep $duration
mv jobIdsMultinodeFull.txt jobIdsMSWeak_ByteEncoder15_Sampling_23

../start_weak.sh ../kurpicz_weak.sh
sleep $duration
mv jobIdsMultinodeFull.txt jobIdsKurpiczWeak

../start_weak.sh ../hQuicksort_weak.sh
sleep $duration
mv jobIdsMultinodeFull.txt jobIdsHQuicksortWeak
