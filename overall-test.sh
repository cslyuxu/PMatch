#####################!/bin/bash

ulimit -s unlimited

#####################label for slashdot and twitter
labelsize="64"

#labelsize="100"

#####################label for dblp
labelsized="64"

#labelsized="150"

hopsize="4"

diameter="3"

queryNum="10"

pathlength="3"

portion="3"

#flag for Baseline of lgpq semantic: 0: sub-iso/hom       1:ssim
flag="1"



##################### Test

#rm -rf results/*

rm CMakeCache.txt

cmake .

make



#for((i=1;i<=$queryNum;i++));
#do
#./main query/Query-Twitter/T$i-$labelsize graph/twitter$labelsize path/Path-Twitter/T$i-$labelsize index/T-Index1-$labelsize index/T-Index2-$labelsize $hopsize
#done







for((i=11;i<=$queryNum;i++));
do
./main query/Query-Slashdot/S$i-$labelsize graph/slashdot$labelsize path/Path-Slashdot/S$i-$labelsize index/S-Index1-$labelsize index/S-Index2-$labelsize $hopsize $pathlength $portion $flag
done

for((i=11;i<=$queryNum;i++));
do
./main query/Query-DBLP/D$i-$labelsized graph/dblp-un$labelsized path/Path-DBLP/D$i-$labelsized index/D-Index1-$labelsized index/D-Index2-$labelsized $hopsize $pathlength $portion $flag
done

for((i=11;i<=$queryNum;i++));
do
./main query/Query-Twitter/T$i-$labelsize graph/twitter$labelsize path/Path-Twitter/T$i-$labelsize index/T-Index1-$labelsize index/T-Index2-$labelsize $hopsize $pathlength $portion $flag
done



labelsize="100"

labelsized="150"

flag="0"

for((i=3;i<=1;i++));
do
./main query/Query-Slashdot/S$i-$labelsize graph/slashdot$labelsize path/Path-Slashdot/S$i-$labelsize index/S-Index1-$labelsize index/S-Index2-$labelsize $hopsize $pathlength $portion $flag
done

for((i=5;i<=6;i++));
do
./main query/Query-DBLP/D$i-$labelsized graph/dblp-un$labelsized path/Path-DBLP/D$i-$labelsized index/D-Index1-$labelsized index/D-Index2-$labelsized $hopsize $pathlength $portion $flag
done

for((i=2;i<=1;i++));
do
./main query/Query-Twitter/T$i-$labelsize graph/twitter$labelsize path/Path-Twitter/T$i-$labelsize index/T-Index1-$labelsize index/T-Index2-$labelsize $hopsize $pathlength $portion $flag
done
