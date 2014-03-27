#!/bin/sh

#CMAKE_PREFIX_PATH=/opt/trilinos/:$CMAKE_PREFIX_PATH \
CXX=mpicxx \
FC=mpif90 \
cmake \
    -DCMAKE_INSTALL_PREFIX:PATH=/opt/nosh/dev/ \
    -DCMAKE_BUILD_TYPE:STRING=Release \
    -DBUILD_SHARED_LIBS:BOOL=TRUE \
    ../../source
