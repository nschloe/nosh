#!/bin/sh

cmake \
  -D CMAKE_INSTALL_PREFIX:PATH=/opt/trilinos/dev/openmpi/1.4.3/gcc/4.7.1/release/shared \
  -D CMAKE_SHARED_LINKER_FLAGS:STRING="-Wl,--no-undefined" \
  -D Trilinos_ENABLE_DEVELOPMENT_MODE:BOOL=OFF \
  -D CMAKE_BUILD_TYPE:STRING=Release \
  -D BUILD_SHARED_LIBS:BOOL=ON \
  -D TPL_FIND_SHARED_LIBS:BOOL=ON \
  -D TPL_ENABLE_MPI:BOOL=ON \
  -D Trilinos_ENABLE_TEUCHOS_TIME_MONITOR:BOOL=ON \
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
