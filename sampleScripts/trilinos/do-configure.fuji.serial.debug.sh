cmake \
  -D CMAKE_C_COMPILER="clang" \
  -D CMAKE_CXX_COMPILER="clang++" \
  -D CMAKE_INSTALL_PREFIX:PATH=/opt/trilinos/dev/master/clang/3.0/debug \
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
  -D Trilinos_ENABLE_SEACASNemslice:BOOL=ON \
  -D Trilinos_ENABLE_SEACASNemspread:BOOL=ON \
  -D TPL_ENABLE_Netcdf:BOOL=ON \
  ../../source
#      -D Netcdf_LIBRARY_DIRS:PATH=/usr/local/ \
#      -D TPL_Netcdf_INCLUDE_DIRS:PATH=/usr/local/ \
