#!/bin/bash
#$1 for hop length

echo "$1 for hop, $2 for label size of slashdot and twitter, $3 for label size of dblp"

rm CMakeCache.txt

cmake .

make

./main graph/slashdot$2 $1 $2

./main graph/dblp-un$3 $1 $3

./main graph/twitter$2 $1 $2

echo "###########################################################################"

echo "###################  Graph Index Construction Done! #######################"

echo "###########################################################################"
