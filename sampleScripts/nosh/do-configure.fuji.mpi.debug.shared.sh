#!/bin/sh

# CMAKE_PREFIX_PATH=/opt/vtk/upstream/:$CMAKE_PREFIX_PATH \
#CMAKE_PREFIX_PATH=/opt/trilinos/upstream/:$CMAKE_PREFIX_PATH \
cmake \
    -DCMAKE_BUILD_TYPE:STRING=Debug \
    -DTIME_MONITOR:BOOL=ON \
    -DCMAKE_CXX_COMPILER:STRING=mpicxx \
    -DBUILD_SHARED_LIBS:BOOL=ON \
    -DCMAKE_SHARED_LINKER_FLAGS:STRING="-Wl,--no-undefined" \
    ../../source/
#    -DCMAKE_INSTALL_PREFIX:PATH=/opt/nosh/dev/openmpi/1.4.5/gcc/4.7.3/debug/shared/ \
