#!/bin/bash

module load gcc/14.2.0
module load python
module load gpu
module load hdf5
module load openmpi
module load cmake

cmake -B/home/lwatan/data/SPH-EXA-fork/build -S/home/lwatan/data/SPH-EXA-fork -DSPH_EXA_WITH_GRACKLE=OFF
cd build/sphexa 
make sphexa