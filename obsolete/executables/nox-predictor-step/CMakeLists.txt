INCLUDE( ../../cmake-system.txt REQUIRED )

ENABLE_LANGUAGE( CXX )

IF(MPI_COMPILED)
    FIND_PACKAGE( MPI REQUIRED )
    SET( LINK_LIBRARIES ${LINK_LIBRARIES}
                        ${MPI_LIBRARY}    )
    INCLUDE_DIRECTORIES ( ${MPI_INCLUDE_PATH} )
ENDIF(MPI_COMPILED)

FIND_PACKAGE( LOCAEpetra REQUIRED )
FIND_PACKAGE( NOXEpetra  REQUIRED )

INCLUDE_DIRECTORIES ( ${Trilinos_INCLUDE_DIR} )

INCLUDE_DIRECTORIES ( ${ginzburg-landau_SOURCE_DIR}/packages/libgl/src
                      ${ginzburg-landau_SOURCE_DIR}/packages/libcartesian/src
                      ${ginzburg-landau_SOURCE_DIR}/packages/libio/src )

SET ( ginzburg_landau_SRCS
      main.cpp )

ADD_EXECUTABLE ( nox-predictor-step.out
                 ${ginzburg_landau_SRCS} )

SET ( internal_LIBS
      gl
      cartesian
      io
    )

# ==============================================================================
# Find other external libs
# -- Almost all of this should actually be handled by Trilinos itself.
FIND_PACKAGE( Boost COMPONENTS filesystem REQUIRED )
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIRS} )
SET( LINK_LIBRARIES ${LINK_LIBRARIES}
                    ${Boost_LIBRARIES} )
# ------------------------------------------------------------------------------
FIND_LIBRARY( GFORTRAN ${GFORTRAN_LIBNAME}
              HINTS ${GFORTRAN_DIR}    )
SET( LINK_LIBRARIES ${LINK_LIBRARIES}
                    ${GFORTRAN} )
# ==============================================================================


# ==============================================================================
# Set executable linking information.
TARGET_LINK_LIBRARIES ( nox-predictor-step.out
                        ${internal_LIBS}
			LOCAEpetra
                        NOXEpetra
                        ${LINK_LIBRARIES}
                      )
# ==============================================================================
