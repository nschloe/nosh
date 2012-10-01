#!/bin/sh

CMAKE_PREFIX_PATH=/opt/trilinos/dev/openmpi/1.4.3/gcc/4.6.3/release/shared:$CMAKE_PREFIX_PATH \
cmake \
    -DCMAKE_CXX_COMPILER=mpicxx \
    -DCMAKE_Fortran_COMPILER=mpif90 \
    -DCMAKE_BUILD_TYPE=Release \
    ../../source
