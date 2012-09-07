#!/bin/bash

module load cmake
module load boost
module load binutils
# Serial netCDF is indeed sufficient.
module load netcdf

#      -D BLAS_LIBRARY_NAMES:STRING=acml \
#      -D BLAS_LIBRARY_DIRS:FILEPATH="$ACML_DIR"/gfortran64/lib/ \
#      -D LAPACK_LIBRARY_NAMES:STRING=acml \
#      -D LAPACK_LIBRARY_DIRS:FILEPATH="$ACML_DIR"/gfortran64/lib/ \

#  -D Trilinos_SKIP_FORTRANCINTERFACE_VERIFY_TEST:BOOL=ON \
cmake \
  -D CMAKE_INSTALL_PREFIX:PATH="$SCRATCH/trilinos/dev/pgi/" \
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
      -D BLAS_LIBRARY_DIRS:PATH="$LIBSCI_BASE_DIR/pgi/119/mc12/lib" \
      -D BLAS_LIBRARY_NAMES="sci_pgi" \
  -D TPL_ENABLE_LAPACK:BOOL=ON \
      -D LAPACK_LIBRARY_DIRS:PATH="$LIBSCI_BASE_DIR/pgi/119/mc12/lib" \
      -D LAPACK_LIBRARY_NAMES="sci_pgi" \
  -D Trilinos_ENABLE_NOX:BOOL=ON \
      -D NOX_ENABLE_LOCA:BOOL=ON \
  -D Trilinos_ENABLE_Piro:BOOL=ON \
  -D Trilinos_ENABLE_ML:BOOL=ON \
  -D Trilinos_ENABLE_Anasazi:BOOL=ON \
  -D Trilinos_ENABLE_STK:BOOL=ON \
  -D Trilinos_ENABLE_Belos:BOOL=ON \
  -D Trilinos_ENABLE_SEACASIoss:BOOL=ON \
  -D TPL_ENABLE_Netcdf:BOOL=ON \
      -D Netcdf_LIBRARY_DIRS:PATH="$CRAY_NETCDF_DIR/pgi/119/lib" \
      -D Netcdf_INCLUDE_DIRS:PATH="$CRAY_NETCDF_DIR/pgi/119/include" \
  ../../source/
