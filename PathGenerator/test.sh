#!/bin/bash
#$1 for diameter $2 for query number

echo "$1 for pathlength, $2 for pathNum, $3 for queryNum, $4 for label size of slashdot and twitter, $5 for label size of dblp"

rm CMakeCache.txt

cmake .

make

for((i=1;i<=$3;i++));
do
./main Query-Slashdot/S$i-$4 $1 $2 $i $4
done


for((i=1;i<=$3;i++));
do
./main Query-DBLP/D$i-$5 $1 $2 $i $5
done

for((i=1;i<=$3;i++));
do
./main Query-Twitter/T$i-$4 $1 $2 $i $4
done
