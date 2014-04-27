#!/bin/bash

command -v cmake >/dev/null 2>&1 || { echo "cmake not found. Try 'module load cmake'." >&2; exit 1; }

# Check which compiler is loadad.
if [ "$CRAY_PRGENVPGI" == "loaded" ]; then
  echo "Using PGI compiler."
elif [ "$CRAY_PRGENVCRAY" == "loaded" ]; then
  echo "Using Cray compiler."
elif [ "$CRAY_PRGENVGNU" == "loaded" ]; then
  echo "Using GCC compiler."
elif [ "$CRAY_PRGENVINTEL" == "loaded" ]; then
  echo "Using Intel compiler."
else
  echo "Unknown compiler suite selected. Abort."
  exit 1
fi

echo ${BOOST_DIR?Error: Need Boost loaded (module load boost).}
#echo ${CRAY_TRILINOS_PREFIX_DIR?Error: Need Trilinos loaded (module load cray-trilinos).}

sleep 5

CMAKE_PREFIX_PATH=$SCRATCH/trilinos-sandybridge:$CMAKE_PREFIX_PATH \
CC=cc \
CXX=CC \
FC=ftn \
cmake \
  -DCMAKE_INSTALL_PREFIX:PATH="$SCRATCH/nosh/dev/intel/" \
  -DCMAKE_BUILD_TYPE:STRING=Release \
  $HOME/software/nosh/dev/
