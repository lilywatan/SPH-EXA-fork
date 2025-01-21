#!/bin/bash
module load gcc/13.2.0
module load python
module load openmpi/5.0.5
module load hdf5/1.12.2
module load cmake

rm -rf build-release
rm -rf build-release/CMakeCache.txt
mkdir build-release
gcc --version
ldd --version 

#export CC=mpicc
#export CXX=mpiCC

cmake -B./build-release -S./ -DSPH_EXA_WITH_GRACKLE=OFF -DCMAKE_BUILD_TYPE=Release

pwd
cd build-release/main/src/sphexa
make clean
make sphexa
cd ../