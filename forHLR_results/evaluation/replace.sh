sed -i -- 's/MS_LCP_S/msLCPS/g'    $1
sed -i -- 's/MS_NoLCP_C/msNoLCPC/g' $1
sed -i -- 's/MS_LCP_C/msLCPC/g' $1
sed -i -- 's/MS_NoLCP_C/msNoLCPC/g' $1
sed -i -- 's/PD_NoGolomb/pdNoGolomb/g' $1
sed -i -- 's/PD_SeqGolomb/pdGolomb/g' $1
sed -i -- 's/PD_NoGolomb/pdNoGolomb/g' $1
sed -i -- 's/PD_SeqGolomb/pdGolomb/g' $1
sed -i -- 's/MS_Simple_S/msSimpleS/g' $1

