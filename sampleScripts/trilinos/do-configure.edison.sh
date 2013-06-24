#!/bin/bash

# Make sure CMake is loaded.
command -v cmake >/dev/null 2>&1 || { echo "cmake not found. Try 'module load cmake'." >&2; exit 1; }

echo ${CRAY_NETCDF_HDF5PARALLEL_VERSION?Error: Need parallel netCDF loaded (module load cray-netcdf-hdf5parallel).}

## Don't use Boost with static libs; this causes unresolved dependencies
## due to errors in Trilinos' build system.
echo ${BOOST_DIR?Error: Need Boost loaded (module load boost).}

# Check which compiler is loadad.
if [ "$CRAY_PRGENVPGI" == "loaded" ]; then
  COMPILER_NAME="pgi"
  SCILIB_NAME="sci_pgi"
  SKIP_CF="ON"
elif [ "$CRAY_PRGENVCRAY" == "loaded" ]; then
  COMPILER_NAME="cray"
  SCILIB_NAME="sci_cray"
  SKIP_CF="OFF"
elif [ "$CRAY_PRGENVGNU" == "loaded" ]; then
  COMPILER_NAME="gnu"
  SCILIB_NAME="sci_gnu"
  SKIP_CF="OFF"
elif [ "$CRAY_PRGENVINTEL" == "loaded" ]; then
  COMPILER_NAME="intel"
  SCILIB_NAME="sci_intel"
  SKIP_CF="OFF"
else
  echo "Unknown compiler suite selected. Abort."
  exit 1
fi

# Give the user some time to ctrl-c out. :)
echo "Using <$COMPILER_NAME> compilers."
sleep 5


# Use TPL_FIND_SHARED_LIBS:BOOL=OFF to force Trilinos to find
# the static versions of the TPLs instead of their shared
# counterparts.

# Do use *static* linking, too.
# Otherwise, one gets errors of the kind
# "/usr/bin/ld: attempted static link of dynamic object".

# Use OpenMP to avoid linking errors of the kind
#
# /opt/cray/libsci/12.0.02/INTEL/130/sandybridge/lib/libsci_intel.a(crayblas_zgemm.o): In function `crayblas_polyzgemm_':
# ../common/src/crayblas_gemm.c:(.text+0x16b): undefined reference to `omp_get_max_threads'
# [...]
#

# On edison, BLAS and LAPACK are best provided by MKL. To use it, add
# -mkl=cluster to the link line.
# This works out of the box for the Intel compilers; all others must
# load the mkl module first.
#
#echo ${CRAY_LIBSCI_PREFIX_DIR?Error: Need netCDF loaded (module load cray-libsci).}

#      -D BLAS_LIBRARY_DIRS:PATH="$CRAY_LIBSCI_PREFIX_DIR/lib" \
#      -D BLAS_LIBRARY_NAMES:STRING="${SCILIB_NAME}" \
#
#      -D LAPACK_LIBRARY_DIRS:PATH="$CRAY_LIBSCI_PREFIX_DIR/lib" \
#      -D LAPACK_LIBRARY_NAMES:STRING="${SCILIB_NAME}" \

cmake \
  -D CMAKE_INSTALL_PREFIX:PATH="$SCRATCH/trilinos/dev/${COMPILER_NAME}/" \
  -D CMAKE_BUILD_TYPE:STRING=Release \
  -D Trilinos_ENABLE_DEVELOPMENT_MODE:BOOL=OFF \
  -D Trilinos_ENABLE_TESTS:BOOL=ON \
  -D BUILD_SHARED_LIBS:BOOL=OFF \
  -D TPL_FIND_SHARED_LIBS:BOOL=OFF \
  -D Trilinos_LINK_SEARCH_START_STATIC:BOOL=ON \
  -D Trilinos_ENABLE_OpenMP:BOOL=ON \
  -D TPL_ENABLE_MPI:BOOL=ON \
      -D MPI_C_COMPILER:FILEPATH=${CRAYPE_DIR}/bin/cc \
      -D MPI_CXX_COMPILER:FILEPATH=${CRAYPE_DIR}/bin/CC \
      -D MPI_Fortran_COMPILER:FILEPATH=${CRAYPE_DIR}/bin/ftn \
  -D Trilinos_SKIP_FORTRANCINTERFACE_VERIFY_TEST:BOOL=${SKIP_CF} \
  -D Trilinos_ENABLE_TEUCHOS_TIME_MONITOR:BOOL=ON \
  -D TPL_ENABLE_BinUtils:BOOL=OFF \
  -D TPL_ENABLE_Boost:BOOL=ON \
      -D Boost_LIBRARY_DIRS:PATH="$BOOST_DIR/libs" \
      -D Boost_INCLUDE_DIRS:PATH="$BOOST_DIR/include" \
  -D TPL_ENABLE_BLAS:BOOL=ON \
      -D TPL_BLAS_LIBRARIES:STRING="-mkl=cluster" \
  -D TPL_ENABLE_LAPACK:BOOL=ON \
      -D TPL_LAPACK_LIBRARIES:STRING="-mkl=cluster" \
  -D Trilinos_ENABLE_NOX:BOOL=ON \
      -D NOX_ENABLE_LOCA:BOOL=ON \
  -D Trilinos_ENABLE_Piro:BOOL=ON \
  -D Trilinos_ENABLE_Belos:BOOL=ON \
  -D Trilinos_ENABLE_ML:BOOL=ON \
  -D Trilinos_ENABLE_Anasazi:BOOL=ON \
  -D Trilinos_ENABLE_STK:BOOL=ON \
    -D Trilinos_ENABLE_SEACASIoss:BOOL=ON \
  -D TPL_ENABLE_Netcdf:BOOL=ON \
      -D Netcdf_INCLUDE_DIRS:PATH="$NETCDF_DIR/include" \
      -D Netcdf_LIBRARY_DIRS:PATH="$NETCDF_DIR/lib" \
  ${HOME}/software/trilinos/dev/source/
