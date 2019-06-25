duration=1
dir="DToNSkewedWeak_$1"
mkdir $dir
cd $dir
#../start_weak.sh ../prefix_doubling_weak_skewed.sh
#sleep $duration
#mv jobIdsMultinodeFull.txt jobIdspd
#
#../start_weak.sh ../merge_sort_weak_skewed.sh
#sleep $duration
#mv jobIdsMultinodeFull.txt jobIdsms
#
#../start_weak.sh ../kurpicz_weak_skewed.sh
#sleep $duration
#mv jobIdsMultinodeFull.txt jobIdskurpicz

../start_weak.sh ../hQuicksort_weak.sh
sleep $duration
mv jobIdsMultinodeFull.txt jobIdshquick
