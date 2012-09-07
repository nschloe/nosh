cmake \
  -D CMAKE_C_COMPILER="gcc" \
  -D CMAKE_CXX_COMPILER="g++" \
  -D CMAKE_INSTALL_PREFIX:PATH=/opt/trilinos/dev/gcc/4.7.1/debug \
  -D CMAKE_BUILD_TYPE=Debug \
  -D BUILD_SHARED_LIBS:BOOL=OFF \
  -D TPL_ENABLE_MPI:BOOL=OFF \
  -D Trilinos_ENABLE_TEUCHOS_TIME_MONITOR=ON \
  -D Trilinos_ENABLE_NOX:BOOL=ON \
      -D NOX_ENABLE_LOCA:BOOL=ON \
  -D Trilinos_ENABLE_Piro:BOOL=ON \
  -D Trilinos_ENABLE_ML:BOOL=ON \
  -D Trilinos_ENABLE_Anasazi:BOOL=ON \
  -D Trilinos_ENABLE_STK:BOOL=ON \
  -D Trilinos_ENABLE_SEACASIoss:BOOL=ON \
  -D Trilinos_ENABLE_SEACASNemslice:BOOL=ON \
  -D Trilinos_ENABLE_SEACASNemspread:BOOL=ON \
  -D TPL_ENABLE_Netcdf:BOOL=ON \
  ../../source
