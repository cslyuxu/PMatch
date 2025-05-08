#!/bin/bash
#$1 for diameter $2 for query number

echo "$1 for diameter, $2 for query number, $3 for label size of slashdot and twitter, $4 for label size of dblp"

rm CMakeCache.txt

cmake .

make

./main graph/slashdot$3 $1 $2 $3

./main graph/dblp-un$4 $1 $2 $4

./main graph/twitter$3 $1 $2 $3

echo "###########################################################################"

echo "######################  Query Generation Done!  ###########################"

echo "###########################################################################"
