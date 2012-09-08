#!/bin/bash

# Make sure the correct compiler is currently loaded.
if [ "$CRAY_PRGENVCRAY" != "loaded" ]; then
echo "Incorrect compiler. Abort."
exit 1
fi

module load cmake/2.8.9
module load boost
module load binutils
module load netcdf

#  -D Trilinos_SKIP_FORTRANCINTERFACE_VERIFY_TEST:BOOL=ON \
cmake \
  -D CMAKE_INSTALL_PREFIX:PATH="$SCRATCH/trilinos/dev/cray/" \
  -D CMAKE_BUILD_TYPE=Release \
  -D Trilinos_ENABLE_DEVELOPMENT_MODE:BOOL=OFF \
  -D TPL_ENABLE_MPI:BOOL=ON \
      -D MPI_C_COMPILER:FILEPATH="$ASYNCPE_DIR/bin/cc" \
      -D MPI_CXX_COMPILER:FILEPATH="$ASYNCPE_DIR/bin/CC" \
      -D MPI_Fortran_COMPILER:FILEPATH="$ASYNCPE_DIR/bin/ftn" \
  -D Trilinos_ENABLE_TEUCHOS_TIME_MONITOR:BOOL=ON \
  -D TPL_ENABLE_Boost:BOOL=ON \
      -D Boost_LIBRARY_DIRS:PATH="$BOOST_DIR/libs" \
      -D Boost_INCLUDE_DIRS:PATH="$BOOST_DIR/include" \
  -D TPL_ENABLE_BLAS:BOOL=ON \
      -D BLAS_LIBRARY_DIRS:PATH="$LIBSCI_BASE_DIR/cray/80/mc12/lib" \
      -D BLAS_LIBRARY_NAMES="sci_cray" \
  -D TPL_ENABLE_LAPACK:BOOL=ON \
      -D LAPACK_LIBRARY_DIRS:PATH="$LIBSCI_BASE_DIR/cray/80/mc12/lib" \
      -D LAPACK_LIBRARY_NAMES="sci_cray" \
  -D Trilinos_ENABLE_NOX:BOOL=ON \
      -D NOX_ENABLE_LOCA:BOOL=ON \
  -D Trilinos_ENABLE_Piro:BOOL=ON \
  -D Trilinos_ENABLE_ML:BOOL=ON \
  -D Trilinos_ENABLE_Anasazi:BOOL=ON \
  -D Trilinos_ENABLE_STK:BOOL=ON \
  -D Trilinos_ENABLE_Belos:BOOL=ON \
  -D Trilinos_ENABLE_SEACASIoss:BOOL=ON \
  -D TPL_ENABLE_Netcdf:BOOL=ON \
      -D Netcdf_LIBRARY_DIRS:PATH="$CRAY_NETCDF_DIR/cray/80/lib" \
      -D Netcdf_INCLUDE_DIRS:PATH="$CRAY_NETCDF_DIR/cray/80/include" \
  ../../source/
