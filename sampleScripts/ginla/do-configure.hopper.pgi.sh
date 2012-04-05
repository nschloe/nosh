#!/bin/bash

#module purge
#module load cmake
#module load boost

#      -D Boost_LIBRARY_DIRS:PATH="$BOOST_DIR"/libs \
#      -D Boost_INCLUDE_DIRS:PATH="$BOOST_DIR" \

#    -D BOOST_INCLUDEDIR:PATH="${BOOST_DIR}" \
#    -D Boost_LIBRARYDIR:PATH="${BOOST_DIR}/boost" \

export CMAKE_PREFIX_PATH="${SCRATCH}/trilinos/10.10.1/pgi/"
cmake \
    -DCMAKE_BUILD_TYPE=Release \
    ../../source/ginla
