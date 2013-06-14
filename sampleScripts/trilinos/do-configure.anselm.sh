#
# configure script for anselm.it4i.cz
#

# Make sure to:
# module load cmake
# module load PrgEnv-intel
# module load intel # for MKL
# module load netcdf

# Make sure CMake is loaded.
command -v cmake >/dev/null 2>&1 || { echo "cmake not found. Try 'module load cmake'." >&2; exit 1; }

if [ -z $MKLROOT ]; then
  echo "MKL not found. Make sure to `module load intel`."
  exit 1;
fi

if [ -z $NETCDF_LIB_DIR ]; then
  echo "NetCDF not found. Make sure to `module load netcdf`."
  exit 1;
fi

# Don't include Sundance for bug <https://software.sandia.gov/bugzilla/show_bug.cgi?id=5927>.

# Don't include threading, cf. <http://comments.gmane.org/gmane.comp.mathematics.libmesh.user/4555>.
#  -D Trilinos_ENABLE_OpenMP:BOOL=ON \

cmake \
  -D CMAKE_INSTALL_PREFIX:PATH=$HOME/trilinos/11.2.3/ \
  -D CMAKE_CXX_FLAGS:STRING="-DMPICH_IGNORE_CXX_SEEK" \
  -D CMAKE_SHARED_LINKER_FLAGS:STRING="-Wl,--no-undefined" \
  -D Trilinos_ENABLE_DEVELOPMENT_MODE:BOOL=OFF \
  -D CMAKE_BUILD_TYPE:STRING=Release \
  -D CMAKE_BUILD_TYPE:STRING="None" \
  -D Trilinos_ENABLE_TESTS:BOOL=ON \
    -D Anasazi_ENABLE_TESTS:BOOL=OFF \
    -D Belos_ENABLE_TESTS:BOOL=OFF \
    -D Sacado_ENABLE_TESTS:BOOL=OFF \
  -D Trilinos_ENABLE_EXAMPLES:BOOL=ON \
  -D BUILD_SHARED_LIBS:BOOL=ON \
  -D TPL_FIND_SHARED_LIBS:BOOL=ON \
  -D TPL_ENABLE_MPI:BOOL=ON \
  -D TPL_ENABLE_BinUtils:BOOL=ON \
  -D TPL_ENABLE_MKL:BOOL=ON \
    -D MKL_LIBRARY_DIRS:FILEPATH="${MKLROOT}/lib/intel64" \
    -D MKL_LIBRARY_NAMES:STRING="mkl_intel_lp64;mkl_sequential;mkl_core" \
    -D MKL_INCLUDE_DIRS:FILEPATH="${MKLROOT}/include" \
  -D TPL_ENABLE_BLAS:BOOL=ON \
    -D BLAS_LIBRARY_DIRS:PATH="${MKLROOT}/lib/intel64" \
    -D BLAS_LIBRARY_NAMES:STRING="mkl_intel_lp64;mkl_sequential;mkl_core" \
  -D TPL_ENABLE_LAPACK:BOOL=ON \
    -D LAPACK_LIBRARY_DIRS:FILEPATH="${MKLROOT}/lib/intel64" \
    -D LAPACK_LIBRARY_NAMES:STRING="mkl_intel_lp64;mkl_sequential;mkl_core" \
  -D TPL_ENABLE_Boost:BOOL=OFF \
      -D Boost_LIBRARY_DIRS:PATH="$BOOST_LIB" \
      -D Boost_INCLUDE_DIRS:PATH="$BOOST_INC" \
  -D TPL_ENABLE_Netcdf:BOOL=ON \
      -D Netcdf_LIBRARY_DIRS:PATH="$NETCDF_LIB_DIR" \
      -D Netcdf_INCLUDE_DIRS:PATH="$NETCDF_INC_DIR" \
  -D Trilinos_ENABLE_ALL_PACKAGES:BOOL=ON \
  -D Trilinos_ENABLE_SECONDARY_STABLE_CODE:BOOL=ON \
  -D Trilinos_ASSERT_MISSING_PACKAGES:BOOL=OFF \
  -D Trilinos_ENABLE_Optika:BOOL=OFF \
  -D Trilinos_ENABLE_PyTrilinos:BOOL=OFF \
  -D Trilinos_ENABLE_Sundance:BOOL=OFF \
  -D Trilinos_ENABLE_TriKota:BOOL=OFF \
  -D TPL_ENABLE_MATLAB:BOOL=OFF \
  -D TPL_ENABLE_Matio:BOOL=OFF \
  -D SEACASExodus_ENABLE_MPI:BOOL=OFF \
  ../trilinos-11.2.3-Source
