#!/bin/bash

#module load cmake
#module load boost

#      -D Boost_LIBRARY_DIRS:PATH="$BOOST_DIR"/libs \
#      -D Boost_INCLUDE_DIRS:PATH="$BOOST_DIR" \
#    -D BOOST_INCLUDEDIR:PATH="${BOOST_DIR}" \
#    -D Boost_LIBRARYDIR:PATH="${BOOST_DIR}/boost" \

export CMAKE_PREFIX_PATH="${SCRATCH}/trilinos/dev/master/gnu/:${CMAKE_PREFIX_PATH}"
cmake \
    -DCMAKE_BUILD_TYPE=Release \
    ../../source/
