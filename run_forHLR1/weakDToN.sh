duration=3
./start_weak.sh prefix_doubling_weak.sh
sleep $duration
mv jobIdsMultinodeFull.txt jobIdsPD_weak_500000_500

./start_weak.sh hQuicksort_weak.sh
sleep $duration
mv jobIdsMultinodeFull.txt jobIdsHQuicksort_weak_500000_500

./start_weak.sh kurpicz_weak.sh
sleep $duration
mv jobIdsMultinodeFull.txt jobIdsKurpicz_weak_500000_500

./start_weak.sh merge_sort_weak.sh
sleep $duration
mv jobIdsMultinodeFull.txt MergeSort_weak_500000_500
