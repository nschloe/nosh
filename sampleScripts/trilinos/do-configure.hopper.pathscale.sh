#!/bin/bash

# Make sure the correct compiler is currently loaded.
if [ "$CRAY_PRGENVPATHSCALE" != "loaded" ]; then
echo "Incorrect compiler. Abort."
exit 1
fi

module load cmake/2.8.9
module load boost/1.50
module load binutils
module load netcdf

# Building shared libraries fails on Hopper:
# Somehow, the linker CC enforced static linking which leads to errors
# of the type
#
# /usr/bin/ld: attempted static link of dynamic object `libnemesis.so'
#
# when executables are linked.
# A Trilinos bug makes Trilnos prefer shared TPL libraries over their
# static counterparts if both are present, independently of BUILD_SHARED_LIBS
# (bug 5704, cf. https://software.sandia.gov/bugzilla/show_bug.cgi?id=5704).
# As a workaround, the full library names are provided here.

cmake \
  -D CMAKE_INSTALL_PREFIX:PATH="$SCRATCH/trilinos/dev/pathscale/" \
  -D CMAKE_BUILD_TYPE=Release \
  -D Trilinos_ENABLE_DEVELOPMENT_MODE:BOOL=OFF \
  -D BUILD_SHARED_LIBS:BOOL=OFF \
  -D TPL_ENABLE_MPI:BOOL=ON \
      -D MPI_C_COMPILER:FILEPATH="$ASYNCPE_DIR/bin/cc" \
      -D MPI_CXX_COMPILER:FILEPATH="$ASYNCPE_DIR/bin/CC" \
      -D MPI_Fortran_COMPILER:FILEPATH="$ASYNCPE_DIR/bin/ftn" \
  -D Trilinos_SKIP_FORTRANCINTERFACE_VERIFY_TEST:BOOL=ON \
  -D Trilinos_ENABLE_TEUCHOS_TIME_MONITOR:BOOL=ON \
  -D TPL_ENABLE_Boost:BOOL=ON \
      -D Boost_LIBRARY_DIRS:PATH="$BOOST_DIR/libs" \
      -D Boost_INCLUDE_DIRS:PATH="$BOOST_DIR/include" \
  -D TPL_ENABLE_BLAS:BOOL=ON \
      -D BLAS_LIBRARY_DIRS:PATH="$CRAY_LIBSCI_PREFIX_DIR/lib" \
      -D BLAS_LIBRARY_NAMES="libsci_cray.a" \
  -D TPL_ENABLE_LAPACK:BOOL=ON \
      -D LAPACK_LIBRARY_DIRS:PATH="$CRAY_LIBSCI_PREFIX_DIR/lib" \
      -D LAPACK_LIBRARY_NAMES="libsci_cray.a" \
  -D Trilinos_ENABLE_NOX:BOOL=ON \
      -D NOX_ENABLE_LOCA:BOOL=ON \
  -D Trilinos_ENABLE_Piro:BOOL=ON \
  -D Trilinos_ENABLE_ML:BOOL=ON \
  -D Trilinos_ENABLE_Anasazi:BOOL=ON \
  -D Trilinos_ENABLE_STK:BOOL=ON \
  -D Trilinos_ENABLE_Belos:BOOL=ON \
  -D Trilinos_ENABLE_SEACASIoss:BOOL=ON \
  -D TPL_ENABLE_Netcdf:BOOL=ON \
      -D Netcdf_INCLUDE_DIRS:PATH="$NETCDF_DIR/include" \
      -D Netcdf_LIBRARY_DIRS:PATH="$NETCDF_DIR/lib" \
      -D Netcdf_LIBRARY_NAMES="libnetcdf.a" \
  ../../source/
