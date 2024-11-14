# CMake generated Testfile for 
# Source directory: /Users/lilywatanabe/Desktop/eth/thesis/SPH-EXA-fork/ryoanji/test/interface
# Build directory: /Users/lilywatanabe/Desktop/eth/thesis/SPH-EXA-fork/cmake-build-release/ryoanji/test/interface
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test([=[global_upsweep_cpu]=] "/opt/homebrew/bin/mpiexec" "-n" "10" "/Users/lilywatanabe/Desktop/eth/thesis/SPH-EXA-fork/cmake-build-release/ryoanji/test/interface/global_upsweep_cpu")
set_tests_properties([=[global_upsweep_cpu]=] PROPERTIES  _BACKTRACE_TRIPLES "/Users/lilywatanabe/Desktop/eth/thesis/SPH-EXA-fork/ryoanji/cmake/ryoanji_add_test.cmake;39;add_test;/Users/lilywatanabe/Desktop/eth/thesis/SPH-EXA-fork/ryoanji/test/interface/CMakeLists.txt;8;ryoanji_add_test;/Users/lilywatanabe/Desktop/eth/thesis/SPH-EXA-fork/ryoanji/test/interface/CMakeLists.txt;12;addRyoanjiMpiTest;/Users/lilywatanabe/Desktop/eth/thesis/SPH-EXA-fork/ryoanji/test/interface/CMakeLists.txt;0;")
