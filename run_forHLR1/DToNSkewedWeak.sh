duration=3
dir="DToNSkewedWeak_$1"
mkdir $dir
cd $dir
../start_weak.sh ../prefix_doubling_weak_skewed.sh
sleep $duration
mv jobIdsMultinodeFull.txt jobIdsPDSkewedWeak_GolombEncoding_01

../start_weak.sh ../merge_sort_weak_skewed.sh
sleep $duration
mv jobIdsMultinodeFull.txt jobIdsMSSkewedWeak_ByteEncoder15_Sampling_23

../start_weak.sh ../kurpicz_weak_skewed.sh
sleep $duration
mv jobIdsMultinodeFull.txt jobIdsKurpiczSkewedWeak

#../start_weak.sh ../hQuicksort_weak.sh
#sleep $duration
#mv jobIdsMultinodeFull.txt jobIdsHQuicksortWeak
