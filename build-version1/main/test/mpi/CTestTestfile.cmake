# CMake generated Testfile for 
# Source directory: /home/lwatan/data/SPH-EXA-fork/main/test/mpi
# Build directory: /home/lwatan/data/SPH-EXA-fork/build-version1/main/test/mpi
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test([=[HDF5IO]=] "/apps/opt/spack/linux-ubuntu20.04-x86_64/gcc-9.3.0/openmpi-4.1.3-efaqdqjshm5pznecwekvbmt42qbwhinl/bin/mpiexec" "-n" "2" "/home/lwatan/data/SPH-EXA-fork/build-version1/main/test/mpi/hdf5io")
set_tests_properties([=[HDF5IO]=] PROPERTIES  _BACKTRACE_TRIPLES "/home/lwatan/data/SPH-EXA-fork/domain/cmake/cstone_add_test.cmake;39;add_test;/home/lwatan/data/SPH-EXA-fork/domain/test/integration_mpi/CMakeLists.txt;10;cstone_add_test;/home/lwatan/data/SPH-EXA-fork/main/test/mpi/CMakeLists.txt;7;addMpiTest;/home/lwatan/data/SPH-EXA-fork/main/test/mpi/CMakeLists.txt;12;addFrontendMpiTest;/home/lwatan/data/SPH-EXA-fork/main/test/mpi/CMakeLists.txt;0;")
