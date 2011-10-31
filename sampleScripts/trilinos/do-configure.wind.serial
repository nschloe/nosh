EXTRA_ARGS=$@
TRILINOS_HOME=../../source

#      -D BLAS_LIBRARY_DIRS:FILEPATH=/opt/atlas/3.9.32/gcc/4.4.5/lib \
#      -D BLAS_LIBRARY_NAMES:STRING="f77blas;cblas;atlas" \
#      -D LAPACK_LIBRARY_DIRS:FILEPATH=/opt/atlas/3.9.32/gcc/4.4.5/lib \
#      -D LAPACK_LIBRARY_NAMES:STRING=lapack \

#      -D BLAS_LIBRARY_NAMES:STRING=blas \
#      -D BLAS_LIBRARY_DIRS:FILEPATH=/opt/netlib-blas/gcc/4.4.5/lib \
#      -D LAPACK_LIBRARY_DIRS:FILEPATH=/opt/lapack/3.3.0/gcc/4.4.5/lib \
#      -D LAPACK_LIBRARY_NAMES:STRING=lapack \

#  -D TPL_ENABLE_BLAS:BOOL=ON \
#      -D BLAS_LIBRARY_DIRS:FILEPATH=/opt/atlas/3.9.32/gcc/4.4.5/lib \
#      -D BLAS_LIBRARY_NAMES:STRING="f77blas;cblas;atlas" \
#  -D TPL_ENABLE_LAPACK:BOOL=ON \
#      -D LAPACK_LIBRARY_DIRS:FILEPATH=/opt/atlas/3.9.32/gcc/4.4.5/lib \
#      -D LAPACK_LIBRARY_NAMES:STRING=lapack \

cmake \
  -D CMAKE_INSTALL_PREFIX:PATH=/opt/trilinos/dev/master/gcc/4.4.5/ \
  -D Trilinos_ENABLE_TEUCHOS_TIME_MONITOR:BOOL=ON \
  -D TPL_ENABLE_MPI:BOOL=OFF \
      -D BLAS_LIBRARY_NAMES:STRING=goto2 \
      -D BLAS_LIBRARY_DIRS:FILEPATH=/opt/gotoblas2/1.13/gcc/4.4.5/serial/lib \
  -D TPL_ENABLE_LAPACK:BOOL=ON \
      -D LAPACK_LIBRARY_NAMES:STRING=goto2 \
      -D LAPACK_LIBRARY_DIRS:FILEPATH=/opt/gotoblas2/1.13/gcc/4.4.5/serial/lib \
  -D TPL_ENABLE_Netcdf:BOOL=ON \
      -D Netcdf_LIBRARY_DIRS:FILEPATH=/opt/netcdf/4.1.2/lib \
      -D Netcdf_INCLUDE_DIRS:FILEPATH=/opt/netcdf/4.1.2/include \
  -D Trilinos_ENABLE_Teuchos:BOOL=ON \
      -D Teuchos_ENABLE_COMPLEX:BOOL=ON \
  -D Trilinos_ENABLE_Belos:BOOL=ON \
      -D Belos_ENABLE_TESTS:BOOL=OFF \
  -D Trilinos_ENABLE_Tpetra:BOOL=ON \
  -D Trilinos_ENABLE_Epetra:BOOL=ON \
  -D Trilinos_ENABLE_EpetraExt:BOOL=ON \
  -D Trilinos_ENABLE_Komplex:BOOL=ON \
  -D Trilinos_ENABLE_Stratimikos:BOOL=ON \
  -D Trilinos_ENABLE_Piro:BOOL=ON \
  -D Trilinos_ENABLE_Trios:BOOL=ON \
    -D Trios_ENABLE_EXODUS:BOOL=ON \
    -D Trios_ENABLE_NEMESIS:BOOL=ON \
    -D Trios_ENABLE_IOSS:BOOL=ON \
  -D Trilinos_ENABLE_STK:BOOL=ON \
  -D Trilinos_ENABLE_Thyra:BOOL=ON \
  -D TPL_ENABLE_HDF5:BOOL=OFF \
  -D Trilinos_ENABLE_EpetraExt:BOOL=ON \
  -D Trilinos_ENABLE_NOX:BOOL=ON \
      -D NOX_ENABLE_LOCA:BOOL=ON \
      -D NOX_ENABLE_TESTS:BOOL=OFF \
      -D NOX_ENABLE_EXAMPLES:BOOL=OFF \
  -D Trilinos_ENABLE_AztecOO:BOOL=ON \
  -D Trilinos_ENABLE_Amesos:BOOL=ON \
      -D Amesos_ENABLE_SuperLU:BOOL=ON \
          -D SuperLU_LIBRARY_DIRS:FILEPATH=/opt/superlu/4.1/gcc/4.4.5/lib \
          -D SuperLU_INCLUDE_DIRS:FILEPATH=/opt/superlu/4.1/gcc/4.4.5/include \
          -D SuperLU_LIBRARY_NAMES:STRING=superlu_4.1 \
      -D Amesos_ENABLE_SuperLUDist:BOOL=OFF \
          -D SuperLUDist_LIBRARY_DIRS:FILEPATH=/opt/superlu-dist/2.0/openmpi/1.4.1/gcc/4.5.0/lib \
          -D SuperLUDist_INCLUDE_DIRS:FILEPATH=/opt/superlu-dist/2.0/openmpi/1.4.1/gcc/4.5.0/include \
      -D TPL_ENABLE_ParMETIS:BOOL=OFF \
          -D ParMETIS_LIBRARY_DIRS:FILEPATH=/opt/parmetis/3.1.1/openmpi/1.4/gcc/4.5.0/lib \
          -D ParMETIS_INCLUDE_DIRS:FILEPATH=/opt/parmetis/3.1.1/openmpi/1.4/gcc/4.5.0/include \
  -D Trilinos_ENABLE_Ifpack:BOOL=ON \
  -D Trilinos_ENABLE_Tifpack:BOOL=ON \
  -D Trilinos_ENABLE_ML:BOOL=ON \
  -D Trilinos_ENABLE_Anasazi:BOOL=ON \
  -D Trilinos_ENABLE_SEACAS:BOOL=ON \
      -D SEACAS_ENABLE_APPLICATIONS:BOOL=ON \
      -D SEACAS_ENABLE_NEM_SLICE:BOOL=ON \
      -D SEACAS_ENABLE_NEM_SPREAD:BOOL=ON \
  -D DART_TESTING_TIMEOUT:STRING=600 \
  $EXTRA_ARGS \
  ${TRILINOS_HOME}
