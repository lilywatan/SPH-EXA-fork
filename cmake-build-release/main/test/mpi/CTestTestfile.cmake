# CMake generated Testfile for 
# Source directory: /Users/lilywatanabe/Desktop/eth/thesis/SPH-EXA-fork/main/test/mpi
# Build directory: /Users/lilywatanabe/Desktop/eth/thesis/SPH-EXA-fork/cmake-build-release/main/test/mpi
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test([=[HDF5IO]=] "/opt/homebrew/bin/mpiexec" "-n" "2" "/Users/lilywatanabe/Desktop/eth/thesis/SPH-EXA-fork/cmake-build-release/main/test/mpi/hdf5io")
set_tests_properties([=[HDF5IO]=] PROPERTIES  _BACKTRACE_TRIPLES "/Users/lilywatanabe/Desktop/eth/thesis/SPH-EXA-fork/domain/cmake/cstone_add_test.cmake;39;add_test;/Users/lilywatanabe/Desktop/eth/thesis/SPH-EXA-fork/domain/test/integration_mpi/CMakeLists.txt;10;cstone_add_test;/Users/lilywatanabe/Desktop/eth/thesis/SPH-EXA-fork/main/test/mpi/CMakeLists.txt;7;addMpiTest;/Users/lilywatanabe/Desktop/eth/thesis/SPH-EXA-fork/main/test/mpi/CMakeLists.txt;12;addFrontendMpiTest;/Users/lilywatanabe/Desktop/eth/thesis/SPH-EXA-fork/main/test/mpi/CMakeLists.txt;0;")
