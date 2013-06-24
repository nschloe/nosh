#!/bin/bash

# Make sure CMake is loaded.
command -v cmake >/dev/null 2>&1 || { echo "cmake not found. Try 'module load cmake'." >&2; exit 1;}

#echo ${NETCDF_DIR?Error Need netCDF loaded (module load netcdf).}
#echo ${BOOST_LIB?Error Need Boost loaded (module load boost).}

# Check which compiler is loadad.
if [ "$CRAY_PRGENVPGI" == "loaded" ]; then
  COMPILER_NAME="pgi"
elif [ "$CRAY_PRGENVCRAY" == "loaded" ]; then
  COMPILER_NAME="cray"
elif [ "$CRAY_PRGENVGNU" == "loaded" ]; then
  COMPILER_NAME="gnu"
elif [ "$CRAY_PRGENVINTEL" == "loaded" ]; then
  COMPILER_NAME="intel"
elif [ "$CRAY_PRGENVPATHSCALE" == "loaded" ]; then
  COMPILER_NAME="pathscale"
else
  echo "Unknown compiler suite selected. Abort."
  exit 1
fi

# Give the user some time to CTRL-C out.
echo "Using <$COMPILER_NAME> compilers."
sleep 5

#module load boost
#      -D Boost_LIBRARY_DIRS:PATH="$BOOST_DIR"/libs \
#      -D Boost_INCLUDE_DIRS:PATH="$BOOST_DIR" \
#    -D BOOST_INCLUDEDIR:PATH="${BOOST_DIR}" \
#    -D Boost_LIBRARYDIR:PATH="${BOOST_DIR}/boost" \

CMAKE_PREFIX_PATH="${SCRATCH}/trilinos/dev/$COMPILER_NAME/:${CMAKE_PREFIX_PATH}" \
cmake \
    -DCMAKE_CXX_COMPILER:FILEPATH="$ASYNCPE_DIR/bin/CC" \
    -DCMAKE_Fortran_COMPILER:FILEPATH="$ASYNCPE_DIR/bin/ftn" \
    -DCMAKE_INSTALL_PREFIX:PATH="$SCRATCH/nosh/dev/$COMPILER_NAME/" \
    -DCMAKE_BUILD_TYPE=Release \
    $HOME/software/nosh/dev/
