
# Make sure to:
# module load cmake
# module load gnu/4.5.2 openmpi/1.4.4--gnu--4.5.2
# module load blas/2007--gnu--4.5.2 lapack/3.3.1--gnu--4.5.2
# module load hdf5/1.8.7_ser--gnu--4.5.2 netcdf/4.1.3--gnu--4.5.2

# Don't activate:
#     * PyTrilinos:
#       numpy isn't installed.
#     * netCDF:
#       No version for default GCC compiler version available.

# Other notes:
#     * Add gfortran to the list of standard linker flags.
#       This is necessary since libblas.a and liblapack.a reference
#       symbols from there.
#     * Tell SEACAS to use *serial* netCDF.
#

cmake \
  -D CMAKE_INSTALL_PREFIX:PATH=$CINECA_SCRATCH/trilinos/11.2.3/ \
  -D CMAKE_SHARED_LINKER_FLAGS:STRING="-Wl,--no-undefined" \
  -D Trilinos_ENABLE_DEVELOPMENT_MODE:BOOL=OFF \
  -D CMAKE_BUILD_TYPE:STRING=Release \
  -D CMAKE_BUILD_TYPE:STRING="None" \
  -D Trilinos_ENABLE_TESTS:BOOL=ON \
    -D Sacado_ENABLE_TESTS:BOOL=OFF \
  -D Trilinos_ENABLE_EXAMPLES:BOOL=ON \
  -D CMAKE_SHARED_LINKER_FLAGS:STRING:="-lgfortran" \
  -D BUILD_SHARED_LIBS:BOOL=ON \
  -D TPL_FIND_SHARED_LIBS:BOOL=ON \
  -D TPL_ENABLE_MPI:BOOL=ON \
  -D TPL_ENABLE_BLAS:BOOL=ON \
    -D BLAS_LIBRARY_DIRS:PATH=$BLAS_LIB \
  -D TPL_ENABLE_LAPACK:BOOL=ON \
    -D LAPACK_LIBRARY_DIRS:PATH=$LAPACK_LIB \
  -D TPL_ENABLE_Boost:BOOL=OFF \
      -D Boost_LIBRARY_DIRS:PATH="$BOOST_LIB" \
      -D Boost_INCLUDE_DIRS:PATH="$BOOST_INC" \
  -D TPL_ENABLE_Netcdf:BOOL=ON \
      -D Netcdf_LIBRARY_DIRS:PATH="$NETCDF_LIB" \
      -D Netcdf_INCLUDE_DIRS:PATH="$NETCDF_INCLUDE" \
  -D Trilinos_ENABLE_ALL_PACKAGES:BOOL=ON \
  -D Trilinos_ENABLE_SECONDARY_STABLE_CODE:BOOL=ON \
  -D Trilinos_ASSERT_MISSING_PACKAGES:BOOL=OFF \
  -D Trilinos_ENABLE_TriKota:BOOL=OFF \
  -D Trilinos_ENABLE_Optika:BOOL=OFF \
  -D Trilinos_ENABLE_PyTrilinos:BOOL=OFF \
  -D TPL_ENABLE_MATLAB:BOOL=OFF \
  -D TPL_ENABLE_Matio:BOOL=OFF \
  -D SEACASExodus_ENABLE_MPI:BOOL=OFF \
  $HOME/trilinos/trilinos-11.2.3-Source
