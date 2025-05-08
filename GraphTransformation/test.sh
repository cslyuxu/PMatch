#!/bin/bash
#$1 for labelsize

echo "$1 for label size of slashdot and twitter, $2 for label size of dblp"

g++ main.cpp -o main -std=c++17

./main graph/slashdot.txt $1 $2

./main graph/dblp-un.txt $1 $2

./main graph/twitter.txt $1 $2

echo "###########################################################################"

echo "#######################  Graph Transformation Done! #######################"

echo "###########################################################################"
