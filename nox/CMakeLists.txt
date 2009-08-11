ENABLE_LANGUAGE( CXX )

# make use of our own modules
SET(GL_CMAKE_DIR "${ginzburg-landau_SOURCE_DIR}/CMakeModules")
SET(CMAKE_MODULE_PATH "${GL_CMAKE_DIR}" ${CMAKE_MODULE_PATH})

FIND_PACKAGE( Boost )

INCLUDE_DIRECTORIES ( ${ginzburg-landau_SOURCE_DIR}/libgl
                      /usr/local/trilinos/090603a/include/
                      /opt/trilinos-090603a/LINUX_MPI_64BIT/include/
                      /opt/mpich2-1.1/include/
                      /opt/boost/1.39.0/include/boost-1_39/ )

SET ( ginzburg_landau_SRCS
      gl-nox.cpp )


ADD_EXECUTABLE ( gl-nox.out
                 ${ginzburg_landau_SRCS} )

SET ( internal_LIBS
      gl
    )

# ==============================================================================
# Find Trilinos libraries
FIND_LIBRARY( TEUCHOS teuchos
              HINTS /usr/local/trilinos/current/lib
                    /opt/trilinos-090603a/LINUX_MPI_64BIT/lib )
FIND_LIBRARY( EPETRA epetra
              HINTS /usr/local/trilinos/current/lib
                    /opt/trilinos-090603a/LINUX_MPI_64BIT/lib )
FIND_LIBRARY( EPETRAEXT epetraext
              HINTS /usr/local/trilinos/current/lib
                    /opt/trilinos-090603a/LINUX_MPI_64BIT/lib )
FIND_PACKAGE( NOX )
# FIND_LIBRARY( NOX nox
#               HINTS /usr/local/trilinos/current/lib
#                     /opt/trilinos-090603a/LINUX_MPI_64BIT/lib )
FIND_LIBRARY( NOXEPETRA noxepetra
              HINTS /usr/local/trilinos/current/lib
                    /opt/trilinos-090603a/LINUX_MPI_64BIT/lib )
FIND_LIBRARY( AZTECOO aztecoo
              HINTS /usr/local/trilinos/current/lib
                    /opt/trilinos-090603a/LINUX_MPI_64BIT/lib )
FIND_LIBRARY( IFPACK ifpack
              HINTS /usr/local/trilinos/current/lib
                    /opt/trilinos-090603a/LINUX_MPI_64BIT/lib )
FIND_LIBRARY( AMESOS amesos
              HINTS /usr/local/trilinos/current/lib
                    /opt/trilinos-090603a/LINUX_MPI_64BIT/lib )
FIND_LIBRARY( ML ml
              HINTS /usr/local/trilinos/current/lib
                    /opt/trilinos-090603a/LINUX_MPI_64BIT/lib )
FIND_LIBRARY( LOCA loca
              HINTS /usr/local/trilinos/current/lib
                    /opt/trilinos-090603a/LINUX_MPI_64BIT/lib )
SET( TRILINOS_LIBS
     ${NOX}
     ${NOXEPETRA}
     ${LOCA}
     ${AZTECOO}
     ${IFPACK}
     ${AMESOS}
     ${ML}
     ${EPETRAEXT}
     ${EPETRA}
     ${TEUCHOS}
     ${PARACONT} )
# ==============================================================================


# ==============================================================================
# Find other external libs
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
              HINTS /usr/lib/gcc/x86_64-pc-linux-gnu/4.3.3
                    /opt/gcc-4.3.3/lib64 )
FIND_LIBRARY( MPICH mpich
              HINTS /opt/mpich2-1.1/lib/ )
FIND_LIBRARY( PTHREAD pthread )
FIND_LIBRARY( HDF5 hdf5 
              HINTS /opt/hdf5-1.6.9/lib/ )
FIND_LIBRARY( Z z )
# ==============================================================================


# ==============================================================================
# Set executable linking information.
TARGET_LINK_LIBRARIES ( gl-nox.out
			${internal_LIBS}
                        ${TRILINOS_LIBS}
                        ${SUPERLU}
                        ${SUPERLUDIST}
                        ${PARMETIS}
                        ${METIS}
#                        ${MPICH}
                        ${HDF5}
                        ${Z}
                        ${LAPACK}
                        ${BLAS}
#                        ${PTHREAD}
                        ${GFORTRAN}
                      )
# ==============================================================================
