#!/bin/bash

#module swap PrgEnv-pgi PrgEnv-gnu
#  -D BinUtils_LIBRARY_DIRS:PATH="$BINUTILS_DIR/lib" \

module load boost
module load binutils
module load netcdf-hdf5parallel

EXTRA_ARGS=$@
TRILINOS_HOME=../../source

cmake \
  -D CMAKE_INSTALL_PREFIX:PATH="$SCRATCH/trilinos/dev/master/gnu/" \
  -D Trilinos_ENABLE_TEUCHOS_TIME_MONITOR:BOOL=ON \
  -D CMAKE_BUILD_TYPE=Release \
  -D TPL_ENABLE_MPI:BOOL=ON \
      -D MPI_C_COMPILER:FILEPATH="$ASYNCPE_DIR/bin/cc" \
      -D MPI_CXX_COMPILER:FILEPATH="$ASYNCPE_DIR/bin/CC" \
      -D MPI_Fortran_COMPILER:FILEPATH="$ASYNCPE_DIR/bin/ftn" \
  -D TPL_ENABLE_BLAS:BOOL=ON \
      -D BLAS_LIBRARY_DIRS:FILEPATH="$LIBSCI_BASE_DIR/gnu/46/mc12/lib" \
      -D BLAS_LIBRARY_NAMES="sci_gnu" \
  -D TPL_ENABLE_LAPACK:BOOL=ON \
      -D LAPACK_LIBRARY_DIRS:FILEPATH="$LIBSCI_BASE_DIR/gnu/46/mc12/lib" \
      -D LAPACK_LIBRARY_NAMES="sci_gnu" \
  -D TPL_ENABLE_Boost:BOOL=ON \
      -D Boost_LIBRARY_DIRS:PATH="$BOOST_DIR/libs" \
      -D Boost_INCLUDE_DIRS:PATH="$BOOST_DIR/include" \
  -D TPL_ENABLE_Netcdf:BOOL=ON \
      -D Netcdf_LIBRARY_DIRS:FILEPATH="$CRAY_NETCDF_DIR/gnu/46/lib" \
      -D Netcdf_INCLUDE_DIRS:FILEPATH="$CRAY_NETCDF_DIR/gnu/46/include" \
  -D Trilinos_ENABLE_NOX:BOOL=ON \
      -D NOX_ENABLE_LOCA:BOOL=ON \
  -D Trilinos_ENABLE_Piro:BOOL=ON \
  -D Trilinos_ENABLE_ML:BOOL=ON \
  -D Trilinos_ENABLE_Belos:BOOL=ON \
  -D Trilinos_ENABLE_STK:BOOL=ON \
  -D Trilinos_ENABLE_Anasazi:BOOL=ON \
  -D Trilinos_ENABLE_SEACASIoss:BOOL=ON \
  -D DART_TESTING_TIMEOUT:STRING=600 \
  $EXTRA_ARGS \
  ${TRILINOS_HOME}
