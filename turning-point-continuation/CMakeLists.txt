INCLUDE( ../cmake-system.txt REQUIRED )

ENABLE_LANGUAGE( CXX )

IF(MPI_COMPILED)
    FIND_PACKAGE( MPI REQUIRED )
    SET( LINK_LIBRARIES ${LINK_LIBRARIES}
                        ${MPI_LIBRARY}    )
    INCLUDE_DIRECTORIES ( ${MPI_INCLUDE_PATH} )
ENDIF(MPI_COMPILED)

FIND_PACKAGE( LOCAEpetra REQUIRED )

INCLUDE_DIRECTORIES ( ${Trilinos_INCLUDE_DIR} )

INCLUDE_DIRECTORIES ( ${ginzburg-landau_SOURCE_DIR}/libgl
                      ${ginzburg-landau_SOURCE_DIR}/libio )

SET ( ginzburg_landau_SRCS
      glTurningPointContinuation.cpp )

ADD_EXECUTABLE ( glTurningPointContinuation.out
                 ${ginzburg_landau_SRCS} )

SET ( internal_LIBS
      gl
      io
    )

# ==============================================================================
# Find other external libs
FIND_LIBRARY( GFORTRAN ${GFORTRAN_LIBNAME}
              HINTS ${GFORTRAN_DIR}    )
SET( LINK_LIBRARIES ${LINK_LIBRARIES}
                    ${GFORTRAN} )
# ==============================================================================


# ==============================================================================
# Set executable linking information.
TARGET_LINK_LIBRARIES ( glTurningPointContinuation.out
                        ${internal_LIBS}
                        LOCAEpetra
                        ${LINK_LIBRARIES}
                      )
# ==============================================================================
