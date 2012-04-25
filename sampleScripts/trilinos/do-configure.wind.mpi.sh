EXTRA_ARGS=$@
TRILINOS_HOME=../../source

cmake \
  -D CMAKE_INSTALL_PREFIX:PATH="/opt/trilinos/dev/master/openmpi/1.4.1/gcc/4.4.5/release" \
  -D Trilinos_VERBOSE_CONFIGURE:BOOL=OFF \
  -D CMAKE_BUILD_TYPE=Release \
  -D TPL_ENABLE_MPI:BOOL=ON \
      -D MPI_C_COMPILER:FILEPATH="$ASYNCPE_DIR/bin/cc" \
      -D MPI_CXX_COMPILER:FILEPATH="$ASYNCPE_DIR/bin/CC" \
      -D MPI_Fortran_COMPILER:FILEPATH="$ASYNCPE_DIR/bin/ftn" \
  -D Trilinos_ENABLE_TEUCHOS_TIME_MONITOR:BOOL=ON \
  -D Trilinos_SKIP_FORTRANCINTERFACE_VERIFY_TEST:BOOL=OFF \
  -D TPL_ENABLE_Boost:BOOL=ON \
  -D TPL_ENABLE_BLAS:BOOL=ON \
      -D BLAS_LIBRARY_DIRS:PATH="/opt/gotoblas2/1.13/gcc/4.4.5/serial/lib" \
      -D BLAS_LIBRARY_NAMES="goto2" \
  -D TPL_ENABLE_LAPACK:BOOL=ON \
      -D LAPACK_LIBRARY_DIRS:PATH="/opt/gotoblas2/1.13/gcc/4.4.5/serial/lib" \
      -D LAPACK_LIBRARY_NAMES="goto2" \
  -D Trilinos_ENABLE_NOX:BOOL=ON \
      -D NOX_ENABLE_LOCA:BOOL=ON \
  -D Trilinos_ENABLE_Piro:BOOL=ON \
  -D Trilinos_ENABLE_ML:BOOL=ON \
  -D Trilinos_ENABLE_Anasazi:BOOL=ON \
  -D Trilinos_ENABLE_STK:BOOL=ON \
  -D Trilinos_ENABLE_Belos:BOOL=ON \
  -D Trilinos_ENABLE_SEACASIoss:BOOL=ON \
  -D TPL_ENABLE_Netcdf:BOOL=ON \
      -D Netcdf_LIBRARY_DIRS:PATH="/opt/netcdf/4.1.2/lib" \
      -D Netcdf_INCLUDE_DIRS:PATH="/opt/netcdf/4.1.2/include" \
  -D DART_TESTING_TIMEOUT:STRING=600 \
  $EXTRA_ARGS \
  ${TRILINOS_HOME}
