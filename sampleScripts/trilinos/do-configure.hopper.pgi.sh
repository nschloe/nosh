#!/bin/bash

#module swap PrgEnv-pgi PrgEnv-gnu # Fortran/C++ doesn't work with PGI
#module swap xt-asyncpe/5.01 xt-asyncpe/5.04 # to make netcdf work
#module load cmake
module load boost
module load binutils
#module load netcdf-hdf5parallel/4.1.2
module load netcdf-hdf5parallel
#module load acml
#module load superlu

EXTRA_ARGS=$@
TRILINOS_HOME=../../trilinos-10.10.1-Source

#export NETCDF_DIR=/opt/cray/netcdf-hdf5parallel/4.0.1.0/netcdf-hdf5parallel-pgi
#      -D BLAS_LIBRARY_DIRS:FILEPATH="$ACML_DIR"/pgi64_mp/lib/ \
#      -D LAPACK_LIBRARY_DIRS:FILEPATH="$ACML_DIR"/pgi64_mp/lib/ \

#      -D BLAS_LIBRARY_NAMES:STRING=acml \
#      -D BLAS_LIBRARY_DIRS:FILEPATH="$ACML_DIR"/gfortran64/lib/ \
#      -D LAPACK_LIBRARY_NAMES:STRING=acml \
#      -D LAPACK_LIBRARY_DIRS:FILEPATH="$ACML_DIR"/gfortran64/lib/ \

#  -D Trilinos_EXTRA_LINK_FLAGS:STRING="$CRAY_HDF5_DIR/hdf5-parallel-gnu/lib/libhdf5.a" \

#  -D BinUtils_LIBRARY_DIRS:FILEPATH="$GCC_PATH/snos/lib64" \

#  -D Trilinos_ENABLE_Fortran:BOOL=OFF \

cmake \
  -D CMAKE_INSTALL_PREFIX:PATH="$SCRATCH/trilinos/10.10.1/pgi/" \
  -D Trilinos_ENABLE_TEUCHOS_TIME_MONITOR:BOOL=ON \
  -D Trilinos_SKIP_FORTRANCINTERFACE_VERIFY_TEST:BOOL=ON \
  -D TPL_ENABLE_MPI:BOOL=ON \
      -D MPI_C_COMPILER:FILEPATH="$ASYNCPE_DIR/bin/cc" \
      -D MPI_CXX_COMPILER:FILEPATH="$ASYNCPE_DIR/bin/CC" \
      -D MPI_Fortran_COMPILER:FILEPATH="$ASYNCPE_DIR/bin/ftn" \
  -D TPL_ENABLE_BLAS:BOOL=ON \
      -D BLAS_LIBRARY_DIRS:FILEPATH="$LIBSCI_BASE_DIR/pgi/109/mc12/lib" \
      -D BLAS_LIBRARY_NAMES="sci_pgi" \
  -D TPL_ENABLE_LAPACK:BOOL=ON \
      -D LAPACK_LIBRARY_DIRS:FILEPATH="$LIBSCI_BASE_DIR/pgi/109/mc12/lib" \
      -D LAPACK_LIBRARY_NAMES="sci_pgi" \
  -D TPL_ENABLE_Boost:BOOL=ON \
      -D Boost_LIBRARY_DIRS:PATH="$BOOST_DIR/libs" \
      -D Boost_INCLUDE_DIRS:PATH="$BOOST_DIR/include" \
  -D TPL_ENABLE_Netcdf:BOOL=ON \
      -D Netcdf_LIBRARY_DIRS:FILEPATH="$CRAY_NETCDF_DIR/pgi/109/lib" \
      -D Netcdf_INCLUDE_DIRS:FILEPATH="$CRAY_NETCDF_DIR/pgi/109/include" \
  -D BinUtils_LIBRARY_DIRS:PATH="$BINUTILS_DIR/lib" \
  -D Trilinos_ENABLE_Teuchos:BOOL=ON \
  -D Trilinos_ENABLE_Belos:BOOL=ON \
  -D Trilinos_ENABLE_Epetra:BOOL=ON \
  -D Trilinos_ENABLE_EpetraExt:BOOL=ON \
  -D Trilinos_ENABLE_Stratimikos:BOOL=ON \
  -D Trilinos_ENABLE_Piro:BOOL=ON \
  -D Trilinos_ENABLE_STK:BOOL=ON \
  -D Trilinos_ENABLE_NOX:BOOL=ON \
      -D NOX_ENABLE_LOCA:BOOL=ON \
  -D Trilinos_ENABLE_ML:BOOL=ON \
  -D Trilinos_ENABLE_Anasazi:BOOL=ON \
  -D Trilinos_ENABLE_SEACASIoss:BOOL=ON \
  -D DART_TESTING_TIMEOUT:STRING=600 \
  $EXTRA_ARGS \
  ${TRILINOS_HOME}
