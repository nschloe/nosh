#!/bin/sh

export CC="mpicc"
export CXX="mpic++"

cmake \
  -D CMAKE_INSTALL_PREFIX:PATH=/opt/trilinos/dev/openmpi/1.4.5/gcc/4.7.2/debug/shared \
  -D CMAKE_SHARED_LINKER_FLAGS:STRING="-Wl,--no-undefined" \
  -D CMAKE_BUILD_TYPE:STRING=Debug \
  -D CMAKE_C_FLAGS_DEBUG:STRING="-Og -g -ggdb -Wall -pedantic -fbounds-check -Wextra -Wstrict-null-sentinel -Wshadow -Woverloaded-virtual -Wsign-compare -fsanitize=address" \
  -D CMAKE_CXX_FLAGS_DEBUG:STRING="-Og -g -ggdb -Wall -pedantic -fbounds-check -Wextra -Wstrict-null-sentinel -Wshadow -Woverloaded-virtual -Weffc++ -Wsign-compare -ansi -std=c++11 -fsanitize=address" \
  -D BUILD_SHARED_LIBS:BOOL=ON \
  -D TPL_FIND_SHARED_LIBS:BOOL=ON \
  -D TPL_ENABLE_MPI:BOOL=ON \
  -D Trilinos_ENABLE_TEUCHOS_TIME_MONITOR:BOOL=ON \
  -D Trilinos_ENABLE_NOX:BOOL=ON \
      -D NOX_ENABLE_LOCA:BOOL=ON \
  -D Trilinos_ENABLE_Piro:BOOL=ON \
  -D Trilinos_ENABLE_ML:BOOL=ON \
  -D Trilinos_ENABLE_Anasazi:BOOL=ON \
  -D Trilinos_ENABLE_Stokhos:BOOL=ON \
  -D Trilinos_ENABLE_FEI:BOOL=OFF \
  -D Trilinos_ENABLE_PyTrilinos:BOOL=ON \
  -D Trilinos_ENABLE_STK:BOOL=ON \
  -D Trilinos_ENABLE_SEACASIoss:BOOL=ON \
  -D Trilinos_ENABLE_SEACASNemslice:BOOL=ON \
  -D Trilinos_ENABLE_SEACASNemspread:BOOL=ON \
  -D TPL_ENABLE_Netcdf:BOOL=ON \
  ../../source


#  -D Trilinos_INSTALL_INCLUDE_DIR:PATH=include/trilinos \
#  -D CMAKE_BUILD_TYPE=None \
