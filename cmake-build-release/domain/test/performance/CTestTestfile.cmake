# CMake generated Testfile for 
# Source directory: /Users/lilywatanabe/Desktop/eth/thesis/SPH-EXA-fork/domain/test/performance
# Build directory: /Users/lilywatanabe/Desktop/eth/thesis/SPH-EXA-fork/cmake-build-release/domain/test/performance
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test([=[scan_perf]=] "/opt/homebrew/bin/mpiexec" "-n" "1" "/Users/lilywatanabe/Desktop/eth/thesis/SPH-EXA-fork/cmake-build-release/domain/test/performance/scan_perf")
set_tests_properties([=[scan_perf]=] PROPERTIES  _BACKTRACE_TRIPLES "/Users/lilywatanabe/Desktop/eth/thesis/SPH-EXA-fork/domain/cmake/cstone_add_test.cmake;39;add_test;/Users/lilywatanabe/Desktop/eth/thesis/SPH-EXA-fork/domain/test/performance/CMakeLists.txt;8;cstone_add_test;/Users/lilywatanabe/Desktop/eth/thesis/SPH-EXA-fork/domain/test/performance/CMakeLists.txt;0;")
