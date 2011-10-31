#!/bin/bash

module load cmake
module load boost
module load binutils
#module load acml
#module load netcdf-hdf5parallel
#module load hdf5-parallel
#module load superlu

EXTRA_ARGS=$@
TRILINOS_HOME=../../source

#export NETCDF_DIR=/opt/cray/netcdf-hdf5parallel/4.0.1.0/netcdf-hdf5parallel-pgi
#      -D BLAS_LIBRARY_DIRS:FILEPATH="$ACML_DIR"/pgi64_mp/lib/ \
#      -D LAPACK_LIBRARY_DIRS:FILEPATH="$ACML_DIR"/pgi64_mp/lib/ \

#      -D BLAS_LIBRARY_NAMES:STRING=acml \
#      -D BLAS_LIBRARY_DIRS:FILEPATH="$ACML_DIR"/gfortran64/lib/ \
#      -D LAPACK_LIBRARY_NAMES:STRING=acml \
#      -D LAPACK_LIBRARY_DIRS:FILEPATH="$ACML_DIR"/gfortran64/lib/ \

#  -D Trilinos_EXTRA_LINK_FLAGS:STRING="$CRAY_HDF5_DIR/hdf5-parallel-gnu/lib/libhdf5.a" \

#  -D BinUtils_LIBRARY_DIRS:FILEPATH="$GCC_PATH/snos/lib64" \

cmake \
  -D CMAKE_INSTALL_PREFIX:PATH="$SCRATCH/trilinos/dev/master/gnu/" \
  -D Trilinos_ENABLE_TEUCHOS_TIME_MONITOR:BOOL=ON \
  -D TPL_ENABLE_MPI:BOOL=ON \
      -D MPI_C_COMPILER:FILEPATH="$ASYNCPE_DIR/bin/cc" \
      -D MPI_CXX_COMPILER:FILEPATH="$ASYNCPE_DIR/bin/CC" \
      -D MPI_Fortran_COMPILER:FILEPATH="$ASYNCPE_DIR/bin/ftn" \
  -D TPL_ENABLE_BLAS:BOOL=ON \
      -D BLAS_LIBRARY_DIRS:FILEPATH="$LIBSCI_BASE_DIR/gnu/46/istanbul/lib" \
      -D BLAS_LIBRARY_NAMES="sci_gnu" \
  -D TPL_ENABLE_LAPACK:BOOL=ON \
      -D LAPACK_LIBRARY_DIRS:FILEPATH="$LIBSCI_BASE_DIR/gnu/46/istanbul/lib" \
      -D LAPACK_LIBRARY_NAMES="sci_gnu" \
  -D TPL_ENABLE_Boost:BOOL=ON \
      -D Boost_LIBRARY_DIRS:PATH="$BOOST_DIR/libs" \
      -D Boost_INCLUDE_DIRS:PATH="$BOOST_DIR/include" \
  -D TPL_ENABLE_Netcdf:BOOL=OFF \
      -D Netcdf_LIBRARY_DIRS:FILEPATH="$CRAY_NETCDF_DIR/netcdf-hdf5parallel-pgi/lib" \
      -D Netcdf_INCLUDE_DIRS:FILEPATH="$CRAY_NETCDF_DIR/netcdf-hdf5parallel-pgi/include" \
  -D BinUtils_LIBRARY_DIRS:PATH="$BINUTILS_DIR/lib" \
  -D TPL_ENABLE_HDF5:BOOL=ON \
  -D Trilinos_ENABLE_Teuchos:BOOL=ON \
      -D Teuchos_ENABLE_COMPLEX:BOOL=ON \
  -D Trilinos_ENABLE_Belos:BOOL=ON \
      -D Belos_ENABLE_TESTS:BOOL=OFF \
  -D Trilinos_ENABLE_Tpetra:BOOL=ON \
  -D Trilinos_ENABLE_Epetra:BOOL=ON \
  -D Trilinos_ENABLE_EpetraExt:BOOL=ON \
  -D Trilinos_ENABLE_Komplex:BOOL=ON \
  -D Trilinos_ENABLE_Stratimikos:BOOL=ON \
  -D Trilinos_ENABLE_Piro:BOOL=ON \
  -D Trilinos_ENABLE_Trios:BOOL=ON \
    -D Trios_ENABLE_EXODUS:BOOL=ON \
    -D Trios_ENABLE_NEMESIS:BOOL=ON \
    -D Trios_ENABLE_IOSS:BOOL=ON \
  -D Trilinos_ENABLE_STK:BOOL=ON \
  -D Trilinos_ENABLE_Thyra:BOOL=ON \
  -D TPL_ENABLE_HDF5:BOOL=OFF \
  -D Trilinos_ENABLE_EpetraExt:BOOL=ON \
  -D Trilinos_ENABLE_NOX:BOOL=ON \
      -D NOX_ENABLE_LOCA:BOOL=ON \
      -D NOX_ENABLE_TESTS:BOOL=OFF \
      -D NOX_ENABLE_EXAMPLES:BOOL=OFF \
  -D Trilinos_ENABLE_AztecOO:BOOL=ON \
  -D Trilinos_ENABLE_Amesos:BOOL=ON \
      -D Amesos_ENABLE_SuperLU:BOOL=OFF \
      -D Amesos_ENABLE_SuperLUDist:BOOL=OFF \
      -D TPL_ENABLE_ParMETIS:BOOL=OFF \
  -D Trilinos_ENABLE_Ifpack:BOOL=ON \
  -D Trilinos_ENABLE_Tifpack:BOOL=ON \
  -D Trilinos_ENABLE_ML:BOOL=ON \
  -D Trilinos_ENABLE_Anasazi:BOOL=ON \
  -D Trilinos_ENABLE_SEACAS:BOOL=OFF \
  -D Trilinos_ENABLE_SEACASIoss:BOOL=OFF \
      -D SEACAS_ENABLE_APPLICATIONS:BOOL=ON \
      -D SEACAS_ENABLE_CONJOIN:BOOL=ON \
      -D SEACAS_ENABLE_NEM_SLICE:BOOL=ON \
      -D SEACAS_ENABLE_NEM_SPREAD:BOOL=ON \
  -D DART_TESTING_TIMEOUT:STRING=600 \
  $EXTRA_ARGS \
  ${TRILINOS_HOME}
