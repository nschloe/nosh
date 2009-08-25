INCLUDE( ../cmake-system.txt REQUIRED )

ENABLE_LANGUAGE( CXX )

SET(GL_CMAKE_DIR "${ginzburg-landau_SOURCE_DIR}/CMakeModules")
SET(CMAKE_MODULE_PATH "${GL_CMAKE_DIR}" ${CMAKE_MODULE_PATH})

INCLUDE_DIRECTORIES ( ${ginzburg-landau_SOURCE_DIR}/libgl
                      /opt/trilinos-090603a/LINUX_MPI_64BIT/include/ 
                      /opt/boost/1.39.0/include/boost-1_39/ )

SET ( ginzburg_landau_SRCS
      gl-loca.cpp )


ADD_EXECUTABLE ( gl-loca.out
                 ${ginzburg_landau_SRCS} )

SET ( internal_LIBS
      gl
    )

FIND_PACKAGE( Trilinos COMPONENTS NOXEpetra
                                  NOX
                                  AztecOO
                                  IFPACK
                                  Amesos
                                  ML
                                  LOCA
                                  LOCAEpetra
                                  EpetraExt
                                  Epetra
                                  Teuchos )

# ==============================================================================
# Find other external libs
FIND_LIBRARY( PARACONT paracont
              HINTS /home/nschloe/svn/paracont/trunk/lib )
IF(TRILINOS_SUPERLU)
    FIND_LIBRARY( SUPERLU superlu_3.0
                  HINTS /usr/local/superlu/3.0/lib
                        /opt/superlu-3.0/lib )
    SET( LINK_LIBRARIES ${LINK_LIBRARIES}
                        ${SUPERLU_LIBRARY} )
ENDIF(TRILINOS_SUPERLU)
IF(TRILINOS_SUPERLUDIST)
    FIND_LIBRARY( SUPERLUDIST superlu_dist_2.3
                  HINTS /opt/superlu-dist-2.3/lib )
    SET( LINK_LIBRARIES ${LINK_LIBRARIES}
                        ${SUPERLUDIST_LIBRARY} )
ENDIF(TRILINOS_SUPERLUDIST)
IF(TRILINOS_PARMETIS)
    FIND_LIBRARY( PARMETIS parmetis
                  HINTS /opt/parmetis-3.1.1/lib )
    SET( LINK_LIBRARIES ${LINK_LIBRARIES}
                        ${PARMETIS_LIBRARY})
ENDIF(TRILINOS_PARMETIS)
IF(TRILINOS_METIS)
    FIND_LIBRARY( METIS metis
                  HINTS /opt/parmetis-3.1.1/lib )
    SET( LINK_LIBRARIES ${LINK_LIBRARIES}
                        ${METIS_LIBRARY})
ENDIF(TRILINOS_METIS)
FIND_LIBRARY( BLAS goto-blas-gfortran
              HINTS /usr/local/lib
                    /opt/goto-blas-1.26/lib )
FIND_LIBRARY( LAPACK lapack-gfortran
              HINTS /usr/local/lib
                    /opt/lapack-3.2.1/lib )
FIND_LIBRARY( GFORTRAN gfortran
              HINTS /usr/lib/gcc/x86_64-pc-linux-gnu/4.3.3
                    /opt/gcc-4.3.3/lib64 )
FIND_LIBRARY( HDF5 hdf5 
              HINTS /opt/hdf5-1.6.9/lib/ )
FIND_LIBRARY( Z z )
# ==============================================================================


# ==============================================================================
# Set executable linking information.
TARGET_LINK_LIBRARIES ( gl-loca.out
			${internal_LIBS}
                        ${TRILINOS_LIBS}
                        ${LINK_LIBRARIES}
                        ${HDF5}
                        ${Z}
                        ${LAPACK}
                        ${BLAS}
                        ${GFORTRAN}
                      )
# ==============================================================================
