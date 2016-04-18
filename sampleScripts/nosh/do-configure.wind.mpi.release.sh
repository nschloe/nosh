#!/bin/sh

# OMPI_CXX=clang++ \
CMAKE_PREFIX_PATH=/opt/trilinos/launchpad/:$CMAKE_PREFIX_PATH \
CXX=mpicxx \
FC=mpif90 \
cmake \
  -DCMAKE_SHARED_LINKER_FLAGS:STRING="-Wl,--no-undefined" \
  -DCMAKE_INSTALL_PREFIX:PATH=/opt/nosh/dev/ \
  -DCMAKE_BUILD_TYPE:STRING=Release \
  -DBUILD_SHARED_LIBS:BOOL=TRUE \
  ../../source
#  -DENABLE_GCOV:BOOL=TRUE \
