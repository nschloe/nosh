EXTRA_ARGS=$@
TRILINOS_HOME=../../source

cmake \
  -D CMAKE_INSTALL_PREFIX:PATH=/opt/trilinos/dev/master/gcc/4.6.3/debug \
  -D CMAKE_BUILD_TYPE=Debug \
  -D TPL_ENABLE_MPI:BOOL=OFF \
  -D Trilinos_ENABLE_TEUCHOS_TIME_MONITOR=ON \
  -D Trilinos_ENABLE_NOX:BOOL=ON \
      -D NOX_ENABLE_LOCA:BOOL=ON \
  -D Trilinos_ENABLE_Piro:BOOL=ON \
  -D Trilinos_ENABLE_ML:BOOL=ON \
  -D Trilinos_ENABLE_Anasazi:BOOL=ON \
  -D Trilinos_ENABLE_STK:BOOL=ON \
  -D Trilinos_ENABLE_SEACASIoss:BOOL=ON \
  -D TPL_ENABLE_Netcdf:BOOL=ON \
      -D Netcdf_LIBRARY_DIRS:PATH=/opt/netcdf/4.1.2/lib/ \
      -D TPL_Netcdf_INCLUDE_DIRS:PATH=/opt/netcdf/4.1.2/include \
  $EXTRA_ARGS \
  ${TRILINOS_HOME}
