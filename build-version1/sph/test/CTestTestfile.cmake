# CMake generated Testfile for 
# Source directory: /home/lwatan/data/SPH-EXA-fork/sph/test
# Build directory: /home/lwatan/data/SPH-EXA-fork/build-version1/sph/test
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test([=[sph_tests]=] "/home/lwatan/data/SPH-EXA-fork/build-version1/sph/test/sph_tests")
set_tests_properties([=[sph_tests]=] PROPERTIES  _BACKTRACE_TRIPLES "/home/lwatan/data/SPH-EXA-fork/sph/test/CMakeLists.txt;16;add_test;/home/lwatan/data/SPH-EXA-fork/sph/test/CMakeLists.txt;0;")
subdirs("hydro_turb")
