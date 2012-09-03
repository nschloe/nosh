cmake \
  -D CMAKE_INSTALL_PREFIX:PATH=/opt/trilinos/dev/master/openmpi/1.4.3/gcc/4.7.1/release \
  -D CMAKE_BUILD_TYPE=Release \
  -D BUILD_SHARED_LIBS:BOOL=OFF \
  -D TPL_ENABLE_MPI:BOOL=ON \
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

  
#  -D Trilinos_INSTALL_INCLUDE_DIR:PATH=include/trilinos \
#  -D CMAKE_BUILD_TYPE=None \
#  -D Trilinos_ENABLE_DEVELOPMENT_MODE:BOOL=OFF \
#  -D CMAKE_CXX_FLAGS="-Werror=format-security" \
#  -D CMAKE_CXX_FLAGS="-Werror=format-security" \
#      -D Netcdf_LIBRARY_DIRS:PATH=/usr/local/ \
#      -D TPL_Netcdf_INCLUDE_DIRS:PATH=/usr/local/ \
