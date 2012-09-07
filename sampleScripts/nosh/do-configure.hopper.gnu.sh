#!/bin/bash

module load cmake/2.8.9
#module load boost

#      -D Boost_LIBRARY_DIRS:PATH="$BOOST_DIR"/libs \
#      -D Boost_INCLUDE_DIRS:PATH="$BOOST_DIR" \
#    -D BOOST_INCLUDEDIR:PATH="${BOOST_DIR}" \
#    -D Boost_LIBRARYDIR:PATH="${BOOST_DIR}/boost" \

CMAKE_PREFIX_PATH="${SCRATCH}/trilinos/dev/gnu/:${CMAKE_PREFIX_PATH}" \
cmake \
    -DCMAKE_INSTALL_PREFIX:PATH="$SCRATCH/nosh/dev/gnu/" \
    -DCMAKE_BUILD_TYPE=Release \
    ../../source/
