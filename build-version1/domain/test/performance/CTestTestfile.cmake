# CMake generated Testfile for 
# Source directory: /home/lwatan/data/SPH-EXA-fork/domain/test/performance
# Build directory: /home/lwatan/data/SPH-EXA-fork/build-version1/domain/test/performance
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test([=[scan_perf]=] "/apps/opt/spack/linux-ubuntu20.04-x86_64/gcc-9.3.0/openmpi-4.1.3-efaqdqjshm5pznecwekvbmt42qbwhinl/bin/mpiexec" "-n" "1" "/home/lwatan/data/SPH-EXA-fork/build-version1/domain/test/performance/scan_perf")
set_tests_properties([=[scan_perf]=] PROPERTIES  _BACKTRACE_TRIPLES "/home/lwatan/data/SPH-EXA-fork/domain/cmake/cstone_add_test.cmake;39;add_test;/home/lwatan/data/SPH-EXA-fork/domain/test/performance/CMakeLists.txt;8;cstone_add_test;/home/lwatan/data/SPH-EXA-fork/domain/test/performance/CMakeLists.txt;0;")
