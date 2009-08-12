ENABLE_LANGUAGE( CXX Fortran )

# make use of our own modules
SET(GL_CMAKE_DIR "${ginzburg-landau_SOURCE_DIR}/CMakeModules")
SET(CMAKE_MODULE_PATH "${GL_CMAKE_DIR}" ${CMAKE_MODULE_PATH})

# directories for wave
SET(CMAKE_INCLUDE_PATH "/opt/boost/1.39.0/include/boost-1_39/"
                       "/opt/trilinos-090603a/LINUX_MPI_64BIT/include/"
   )

FIND_PACKAGE( Boost REQUIRED )
FIND_PACKAGE( MPI REQUIRED )

# FIND_PACKAGE( LAPACK )
FIND_PACKAGE( Trilinos COMPONENTS NOXEpetra
                                  NOX
                                  AztecOO
                                  IFPACK
                                  Amesos
                                  ML
                                  LOCA
                                  EpetraExt
                                  Epetra
                                  Teuchos )

INCLUDE_DIRECTORIES ( ${ginzburg-landau_SOURCE_DIR}/libgl
                      ${Trilinos_INCLUDE_DIR}
                      ${MPI_INCLUDE_PATH}
                      ${Boost_INCLUDE_DIRS} )

SET ( ginzburg_landau_SRCS
      gl-nox.cpp )

ADD_EXECUTABLE ( gl-nox.out
                 ${ginzburg_landau_SRCS} )

SET ( internal_LIBS
      gl
    )

# ==============================================================================
# Find other external libs
FIND_LIBRARY( PARACONT paracont
              HINTS /home/nschloe/svn/paracont/trunk/lib )
FIND_LIBRARY( SUPERLU superlu_3.0
              HINTS /usr/local/superlu/3.0/lib
                    /opt/superlu-3.0/lib )
FIND_LIBRARY( SUPERLUDIST superlu_dist_2.3
              HINTS /opt/superlu-dist-2.3/lib )
FIND_LIBRARY( PARMETIS parmetis
              HINTS /opt/parmetis-3.1.1/lib )
FIND_LIBRARY( METIS metis
              HINTS /opt/parmetis-3.1.1/lib )
FIND_LIBRARY( BLAS goto-blas-gfortran
              HINTS /usr/local/lib
                    /opt/goto-blas-1.26/lib )
FIND_LIBRARY( LAPACK lapack-gfortran
              HINTS /usr/local/lib
                    /opt/lapack-3.2.1/lib )
FIND_LIBRARY( GFORTRAN gfortran
              HINTS /usr/lib/gcc/x86_64-pc-linux-gnu/4.3.4
                    /opt/gcc-4.3.4/lib64 )
FIND_LIBRARY( PTHREAD pthread )
FIND_LIBRARY( HDF5 hdf5
              HINTS /opt/hdf5-1.6.9/lib/ )
FIND_LIBRARY( Z z )
# ==============================================================================


# ==============================================================================
# Set executable linking information.
TARGET_LINK_LIBRARIES ( gl-nox.out
			${internal_LIBS}
                        ${Trilinos_LIBRARIES}
                        ${PARACONT}
                        ${SUPERLU}
                        ${SUPERLUDIST}
                        ${PARMETIS}
                        ${METIS}
                        ${MPI_LIBRARY}
                        ${HDF5}
                        ${Z}
                        ${LAPACK}
                        ${BLAS}
#                        ${PTHREAD}
                        ${GFORTRAN}
                      )
# ==============================================================================
