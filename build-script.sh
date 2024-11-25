#!/bin/bash

module load gcc/14.2.0
module load python
module load gpu
module load openmpi
module load hdf5
module load cmake

rm -rf build
mkdir build-release

cmake -B/home/lwatan/data/SPH-EXA-fork/build-release -S/home/lwatan/data/SPH-EXA-fork -DSPH_EXA_WITH_GRACKLE=OFF -DCMAKE_BUILD_TYPE=Release
cd build-release/main/src/sphexa
make sphexa
cd /home/lwatan/data/SPH-EXA-fork