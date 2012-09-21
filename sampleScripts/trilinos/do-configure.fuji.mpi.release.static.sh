#!/bin/sh

# From packages/seacas/README:
#
# '''
# If using a netcdf library that is built with --enable-netcdf-4,
# then you also need to tell it to link with the hdf5
# libraries. There is not really a good way to do this; the best
# option is:
#
#    -D Trilinos_EXTRA_LINK_FLAGS="-L{path_to_hdf5_libraries} -lhdf5_hl -lhdf5 -lz -lm"
# '''
#
# Additionally, if netCDF was compiled with DAP support, adding libcurl is necessary
# as well.
# This of course is necessary only for statically linked libraries.

# Okay, okay.
# Ideally, we would like to link against the static version of the TPLs using
#  -D TPL_FIND_SHARED_LIBS:BOOL=OFF.
# However, this is just asking for trouble at this point (2012/09/15).
# One of the problems is with Pthreads (Sandia bug https://software.sandia.gov/bugzilla/show_bug.cgi?id=5715),
# another with netCDF's possible dependency on HDF5,
#  -D Trilinos_EXTRA_LINK_FLAGS:STRING="-lhdf5_hl -lhdf5 -lz -lm -lcurl".
cmake \
  -D CMAKE_INSTALL_PREFIX:PATH=/opt/trilinos/dev/openmpi/1.4.3/gcc/4.7.1/release/static \
  -D Trilinos_ENABLE_DEVELOPMENT_MODE:BOOL=OFF \
  -D CMAKE_BUILD_TYPE:STRING=Release \
  -D BUILD_SHARED_LIBS:BOOL=OFF \
  -D TPL_ENABLE_MPI:BOOL=ON \
  -D Trilinos_ENABLE_TEUCHOS_TIME_MONITOR:BOOL=ON \
  -D Trilinos_ENABLE_NOX:BOOL=ON \
      -D NOX_ENABLE_LOCA:BOOL=ON \
  -D Trilinos_ENABLE_Piro:BOOL=ON \
  -D Trilinos_ENABLE_ML:BOOL=ON \
  -D Trilinos_ENABLE_Anasazi:BOOL=ON \
  -D Trilinos_ENABLE_Stokhos:BOOL=ON \
  -D Trilinos_ENABLE_STK:BOOL=ON \
  -D Trilinos_ENABLE_SEACASIoss:BOOL=ON \
  -D Trilinos_ENABLE_SEACASNemslice:BOOL=ON \
  -D Trilinos_ENABLE_SEACASNemspread:BOOL=ON \
  -D TPL_ENABLE_Netcdf:BOOL=ON \
  ../../source

#  -D Trilinos_INSTALL_INCLUDE_DIR:PATH=include/trilinos \
#  -D CMAKE_BUILD_TYPE=None \
