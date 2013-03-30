#!/bin/sh

CXX=mpicxx \
FC=mpif90 \
CMAKE_PREFIX_PATH=/opt/trilinos/dev/openmpi/1.4.3/gcc/4.6.3/release/shared:$CMAKE_PREFIX_PATH \
cmake \
    -DCMAKE_INSTALL_PREFIX:PATH=/opt/nosh/dev/ \
    -DCMAKE_BUILD_TYPE=Release \
    ../../source
