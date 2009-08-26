INCLUDE( ../cmake-system.txt REQUIRED )

ENABLE_LANGUAGE( CXX )

# make use of our own modules
SET(GL_CMAKE_DIR "${ginzburg-landau_SOURCE_DIR}/CMakeModules")
SET(CMAKE_MODULE_PATH "${GL_CMAKE_DIR}" ${CMAKE_MODULE_PATH})

FIND_PACKAGE( Boost REQUIRED )


SET( LINK_LIBRARIES )

IF(MPI_COMPILED)
    FIND_PACKAGE( MPI REQUIRED )
    SET( LINK_LIBRARIES ${LINK_LIBRARIES}
                        ${MPI_LIBRARY}    )
ENDIF(MPI_COMPILED)

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
                      ${PARACONT_DIR}/include
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
              HINTS ${PARACONT_DIR}/lib )
# ------------------------------------------------------------------------------
IF(TRILINOS_SUPERLU)
    FIND_LIBRARY( SUPERLU ${SUPERLU_LIBNAME}
                  HINTS ${SUPERLU_DIR} )
    SET( LINK_LIBRARIES ${LINK_LIBRARIES}
                        ${SUPERLU_LIBRARY} )
ENDIF(TRILINOS_SUPERLU)
# ------------------------------------------------------------------------------
IF(TRILINOS_SUPERLUDIST)
    FIND_LIBRARY( SUPERLUDIST ${SUPERLUDIST_LIBNAME}
                  HINTS ${SUPERLUDIST_DIR} )
    SET( LINK_LIBRARIES ${LINK_LIBRARIES}
                        ${SUPERLUDIST_LIBRARY} )
ENDIF(TRILINOS_SUPERLUDIST)
# ------------------------------------------------------------------------------
IF(TRILINOS_PARMETIS)
    FIND_LIBRARY( PARMETIS ${PARMETIS_LIBNAME}
                  HINTS ${PARMETIS_DIR} )
    SET( LINK_LIBRARIES ${LINK_LIBRARIES}
                        ${PARMETIS_LIBRARY} )
ENDIF(TRILINOS_PARMETIS)
# ------------------------------------------------------------------------------
IF(TRILINOS_METIS)
    FIND_LIBRARY( METIS ${METIS_LIBNAME}
                  HINTS ${METIS_DIR} )
    SET( LINK_LIBRARIES ${LINK_LIBRARIES}
                        ${METIS_LIBRARY} )
ENDIF(TRILINOS_METIS)
# ------------------------------------------------------------------------------
FIND_LIBRARY( BLAS ${BLAS_LIBNAME}
              HINTS ${BLAS_DIR} )
# ------------------------------------------------------------------------------
FIND_LIBRARY( LAPACK ${LAPACK_LIBNAME}
              HINTS ${LAPACK_DIR} )
# ------------------------------------------------------------------------------
FIND_LIBRARY( GFORTRAN ${GFORTRAN_LIB}
              HINTS ${GFORTRAN_DIR}    )
# ------------------------------------------------------------------------------
# FIND_LIBRARY( PTHREAD pthread )
# ------------------------------------------------------------------------------
FIND_LIBRARY( HDF5 ${HDF5_LIBNAME}
              HINTS ${HDF5_DIR} )
# ------------------------------------------------------------------------------
# FIND_LIBRARY( Z z )
# ==============================================================================


# ==============================================================================
# Set executable linking information.
TARGET_LINK_LIBRARIES ( gl-nox.out
			${internal_LIBS}
                        ${Trilinos_LIBRARIES}
                        ${LINK_LIBRARIES}
                        ${PARACONT}
                        ${HDF5}
#                         ${Z}
                        ${LAPACK}
                        ${BLAS}
#                        ${PTHREAD}
                        ${GFORTRAN}
                      )
# ==============================================================================