duration=3
dir="File_$1"
mkdir $dir
path="../$1"
cd $dir
mkdir tmpScripts

../startjobsFromFileIntern.sh prefix_doubling_file.sh $path $dir
#sleep $duration
#mv jobIdsMultinodeFull.txt jobIdsPDFile_GolombEncoding_01

../startjobsFromFileIntern.sh /merge_sort_file.sh $path $dir
#sleep $duration
#mv jobIdsMultinodeFull.txt jobIdsMSFile_ByteEncoder15_Sampling_23

../startjobsFromFileIntern.sh /kurpicz_file.sh $path $dir
#sleep $duration
#mv jobIdsMultinodeFull.txt jobIdsKurpiczFile

../startjobsFromFileIntern.sh /hQuicksort_file.sh $path $dir
#sleep $duration
#mv jobIdsMultinodeFull.txt jobIdsHQuicksortFile
