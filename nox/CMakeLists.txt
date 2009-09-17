INCLUDE( ../cmake-system.txt REQUIRED )

ENABLE_LANGUAGE( CXX )

IF(MPI_COMPILED)
    FIND_PACKAGE( MPI REQUIRED )
    SET( LINK_LIBRARIES ${LINK_LIBRARIES}
                        ${MPI_LIBRARY}    )
    INCLUDE_DIRECTORIES ( ${MPI_INCLUDE_PATH} )
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
				  Anasazi
                                  Teuchos )
INCLUDE_DIRECTORIES ( ${Trilinos_INCLUDE_DIR} )

INCLUDE_DIRECTORIES ( ${ginzburg-landau_SOURCE_DIR}/libgl )

SET ( ginzburg_landau_SRCS
      gl-nox.cpp )

ADD_EXECUTABLE ( gl-nox.out
                 ${ginzburg_landau_SRCS} )

SET ( internal_LIBS
      gl
    )

# ==============================================================================
# Find other external libs
# -- Almost all of this should actually be handled by Trilinos itself.
FIND_PACKAGE( Boost REQUIRED )
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIRS} )
# ------------------------------------------------------------------------------
FIND_LIBRARY( PARACONT paracont
              HINTS ${PARACONT_DIR}/lib )
SET( LINK_LIBRARIES ${LINK_LIBRARIES}
                    ${PARACONT} )
# ------------------------------------------------------------------------------
IF(TRILINOS_SUPERLU)
    FIND_LIBRARY( SUPERLU ${SUPERLU_LIBNAME}
                  HINTS ${SUPERLU_DIR} )
    SET( LINK_LIBRARIES ${LINK_LIBRARIES}
                        ${SUPERLU} )
ENDIF(TRILINOS_SUPERLU)
# ------------------------------------------------------------------------------
IF(TRILINOS_SUPERLUDIST)
    FIND_LIBRARY( SUPERLUDIST ${SUPERLUDIST_LIBNAME}
                  HINTS ${SUPERLUDIST_DIR} )
    SET( LINK_LIBRARIES ${LINK_LIBRARIES}
                        ${SUPERLUDIST} )
ENDIF(TRILINOS_SUPERLUDIST)
# ------------------------------------------------------------------------------
IF(TRILINOS_PARMETIS)
    FIND_LIBRARY( PARMETIS ${PARMETIS_LIBNAME}
                  HINTS ${PARMETIS_DIR} )
    SET( LINK_LIBRARIES ${LINK_LIBRARIES}
                        ${PARMETIS} )
ENDIF(TRILINOS_PARMETIS)
# ------------------------------------------------------------------------------
IF(TRILINOS_METIS)
    FIND_LIBRARY( METIS ${METIS_LIBNAME}
                  HINTS ${METIS_DIR} )
    SET( LINK_LIBRARIES ${LINK_LIBRARIES}
                        ${METIS} )
ENDIF(TRILINOS_METIS)
# ------------------------------------------------------------------------------
FIND_LIBRARY( HDF5 ${HDF5_LIBNAME}
              HINTS ${HDF5_DIR}/lib )
SET( LINK_LIBRARIES ${LINK_LIBRARIES}
                    ${HDF5} )
# ------------------------------------------------------------------------------
FIND_LIBRARY( Z z )
SET( LINK_LIBRARIES ${LINK_LIBRARIES}
                    ${Z} )
# ------------------------------------------------------------------------------
FIND_LIBRARY( LAPACK ${LAPACK_LIBNAME}
              HINTS ${LAPACK_DIR} )
SET( LINK_LIBRARIES ${LINK_LIBRARIES}
                    ${LAPACK} )
# ------------------------------------------------------------------------------
FIND_LIBRARY( BLAS ${BLAS_LIBNAME}
              HINTS ${BLAS_DIR} )
SET( LINK_LIBRARIES ${LINK_LIBRARIES}
                    ${BLAS} )
# ------------------------------------------------------------------------------
# FIND_LIBRARY( PTHREAD pthread )
# ------------------------------------------------------------------------------
FIND_LIBRARY( GFORTRAN ${GFORTRAN_LIBNAME}
              HINTS ${GFORTRAN_DIR}    )
SET( LINK_LIBRARIES ${LINK_LIBRARIES}
                    ${GFORTRAN} )
# ==============================================================================


# ==============================================================================
# Set executable linking information.
TARGET_LINK_LIBRARIES ( gl-nox.out
			${internal_LIBS}
                        ${Trilinos_LIBRARIES}
                        ${LINK_LIBRARIES}
                      )
# ==============================================================================
