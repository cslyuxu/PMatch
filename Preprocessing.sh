#####################!/bin/bash

#####################label for slashdot and twitter
labelsize="64"

#labelsize="100"

#####################label for dblp
labelsized="64"

#labelsized="150"

##### including index <= 6
hopsize="6"

diameter="3"

queryNum="10"

pathlength="3"

#### can ignore this parameter, have not been used
pathNum="2"

##################### GrapghTransformation

cd GraphTransformation

./test.sh $labelsize $labelsized

cd ../

cp -f GraphTransformation/graph/slashdot$labelsize GraphIndex/graph/slashdot$labelsize
cp -f GraphTransformation/graph/dblp-un$labelsized GraphIndex/graph/dblp-un$labelsized
cp -f GraphTransformation/graph/twitter$labelsize GraphIndex/graph/twitter$labelsize
cp -f GraphTransformation/graph/slashdot$labelsize QueryGenerator/graph/slashdot$labelsize
cp -f GraphTransformation/graph/dblp-un$labelsized QueryGenerator/graph/dblp-un$labelsized
cp -f GraphTransformation/graph/twitter$labelsize QueryGenerator/graph/twitter$labelsize
cp -f GraphTransformation/graph/slashdot$labelsize graph/slashdot$labelsize
cp -f GraphTransformation/graph/dblp-un$labelsized graph/dblp-un$labelsized
cp -f GraphTransformation/graph/twitter$labelsize graph/twitter$labelsize


##################### GraphIndex

cd GraphIndex

./test.sh $hopsize $labelsize $labelsized

cd ../

cp -f GraphIndex/graph/S-Index1-$labelsize index/S-Index1-$labelsize
cp -f GraphIndex/graph/S-Index2-$labelsize index/S-Index2-$labelsize
cp -f GraphIndex/graph/D-Index1-$labelsized index/D-Index1-$labelsized
cp -f GraphIndex/graph/D-Index2-$labelsized index/D-Index2-$labelsized
cp -f GraphIndex/graph/T-Index1-$labelsize index/T-Index1-$labelsize
cp -f GraphIndex/graph/T-Index2-$labelsize index/T-Index2-$labelsize



##################### QueryGenerator

cd QueryGenerator

./test.sh $diameter $queryNum $labelsize $labelsized

cd ../

cp -rf QueryGenerator/Query-Slashdot/* PathGenerator/Query-Slashdot/

cp -rf QueryGenerator/Query-DBLP/* PathGenerator/Query-DBLP/

cp -rf QueryGenerator/Query-Twitter/* PathGenerator/Query-Twitter/

cp -rf QueryGenerator/Query-Slashdot/* query/Query-Slashdot/

cp -rf QueryGenerator/Query-DBLP/* query/Query-DBLP/

cp -rf QueryGenerator/Query-Twitter/* query/Query-Twitter/


##################### PathGenerator

cd PathGenerator

./test.sh $pathlength $pathNum $queryNum $labelsize $labelsized

cd ../

cp -rf PathGenerator/Path-Slashdot/* path/Path-Slashdot/

cp -rf PathGenerator/Path-DBLP/* path/Path-DBLP/

cp -rf PathGenerator/Path-Twitter/* path/Path-Twitter/


