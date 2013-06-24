#!/bin/bash

command -v cmake >/dev/null 2>&1 || { echo "cmake not found. Try 'module load cmake'." >&2; exit 1; }

echo ${BOOST_DIR?Error: Need Boost loaded (module load boost).}

if [ "$CRAY_PRGENVINTEL" != "loaded" ]; then
  echo "Use the Intel compiler!"
  exit 1
fi

#module load boost

#      -D Boost_LIBRARY_DIRS:PATH="$BOOST_DIR"/libs \
#      -D Boost_INCLUDE_DIRS:PATH="$BOOST_DIR" \
#    -D BOOST_INCLUDEDIR:PATH="${BOOST_DIR}" \
#    -D Boost_LIBRARYDIR:PATH="${BOOST_DIR}/boost" \

CC=cc \
CXX=CC \
FC=ftn \
CMAKE_PREFIX_PATH="${SCRATCH}/trilinos/dev/intel/:${CMAKE_PREFIX_PATH}" \
cmake \
    -DCMAKE_INSTALL_PREFIX:PATH="$SCRATCH/nosh/dev/intel/" \
    -DCMAKE_BUILD_TYPE=Release \
    $HOME/software/nosh/dev/
