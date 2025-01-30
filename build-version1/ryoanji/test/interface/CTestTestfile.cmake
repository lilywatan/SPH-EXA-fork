# CMake generated Testfile for 
# Source directory: /home/lwatan/data/SPH-EXA-fork/ryoanji/test/interface
# Build directory: /home/lwatan/data/SPH-EXA-fork/build-version1/ryoanji/test/interface
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test([=[global_upsweep_cpu]=] "/apps/opt/spack/linux-ubuntu20.04-x86_64/gcc-9.3.0/openmpi-4.1.3-efaqdqjshm5pznecwekvbmt42qbwhinl/bin/mpiexec" "-n" "10" "/home/lwatan/data/SPH-EXA-fork/build-version1/ryoanji/test/interface/global_upsweep_cpu")
set_tests_properties([=[global_upsweep_cpu]=] PROPERTIES  _BACKTRACE_TRIPLES "/home/lwatan/data/SPH-EXA-fork/ryoanji/cmake/ryoanji_add_test.cmake;39;add_test;/home/lwatan/data/SPH-EXA-fork/ryoanji/test/interface/CMakeLists.txt;8;ryoanji_add_test;/home/lwatan/data/SPH-EXA-fork/ryoanji/test/interface/CMakeLists.txt;12;addRyoanjiMpiTest;/home/lwatan/data/SPH-EXA-fork/ryoanji/test/interface/CMakeLists.txt;0;")
