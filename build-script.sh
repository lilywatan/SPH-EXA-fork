#!/bin/bash
module load stack/.2024-06-silent 
module load gcc/12.2.0
module load python
module load openmpi/4.1.6
module load hdf5/1.14.3
module load cmake

rm -rf build-release
mkdir build-release
mpicc --version
ldd --version 

export CC=mpicc
export CXX=mpiCC

cmake -B./build-release -S./ -DSPH_EXA_WITH_GRACKLE=OFF -DCMAKE_BUILD_TYPE=Release -DCMAKE_VERBOSE_MAKEFILE=ON 

pwd
cd build-release/main/src/sphexa
make clean
make sphexa
cd ../