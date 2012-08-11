# Example CMake configuration line for an external program that links against Nosh.
# The sources (including CMakeLists.txt) are in ../source (but can of course sit
# anywhere), and the installation path of Nosh is specified in the CMAKE_PREFIX_PATH.
CMAKE_PREFIX_PATH=/opt/nosh/dev/master/openmpi/1.4.3/gcc/4.6.3/release/:$CMAKE_PREFIX_PATH \
cmake ../source
