#!/bin/bash

# Define build directory name (default to 'build-release' if not provided)
BUILD_DIR=${1:-build-release}

echo "Compiling in directory: $BUILD_DIR"

# Load required modules
module load gcc/13.2.0
module load python
module load openmpi/5.0.5
module load hdf5/1.12.2
module load cmake

# Remove previous compilation files safely
rm -rf "$BUILD_DIR"
mkdir "$BUILD_DIR"

gcc --version
ldd --version 

# Uncomment these if using MPI compiler wrappers
# export CC=mpicc
# export CXX=mpiCC

# Configure the project with a separate build folder
cmake -B"$BUILD_DIR" -S./ -DSPH_EXA_WITH_GRACKLE=OFF -DCMAKE_BUILD_TYPE=Release

# Move into the build directory and compile
cd "$BUILD_DIR/main/src/sphexa" || exit
make clean
make sphexa

# Return to the original directory
cd ../../../
echo "Compilation finished. Executable is in $BUILD_DIR/main/src/sphexa/"
