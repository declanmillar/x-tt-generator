# root-tuple

Based on RootTuple https://roottuple.hepforge.org, with an additional TTree for storing process information, e.g cross section, beam info etc.

git clone ssh://git@gitlab.cern.ch:7999/demillar/root-tuple.git

cd root-tuple

mkdir build

cd build

cmake ..

make

cp src/libRootTuple.* <lib/dir>