ENABLE_LANGUAGE( CXX )

IF(MPI_COMPILED)
    FIND_PACKAGE( MPI REQUIRED )
    SET( LINK_LIBRARIES ${LINK_LIBRARIES}
                        ${MPI_LIBRARY}    )
    INCLUDE_DIRECTORIES ( ${MPI_INCLUDE_PATH} )
ENDIF(MPI_COMPILED)

FIND_PACKAGE( LOCA      REQUIRED )
FIND_PACKAGE( NOXEpetra REQUIRED )

FIND_PACKAGE( VTK REQUIRED )
INCLUDE_DIRECTORIES( ${VTK_INCLUDE_DIRS} )
LINK_DIRECTORIES( ${VTK_LIBRARY_DIRS} )

FIND_PACKAGE( Boost REQUIRED )
LINK_DIRECTORIES( ${Boost_LIBRARY_DIRS} )

INCLUDE_DIRECTORIES ( ${Trilinos_INCLUDE_DIR} )

INCLUDE_DIRECTORIES ( ${ginzburg-landau_SOURCE_DIR}/packages/libio/src
                      ${ginzburg-landau_SOURCE_DIR}/packages/libcartesian/src
                      ${ginzburg-landau_SOURCE_DIR}/packages/libgl/src )

SET ( ginzburg_landau_SRCS
      glNoxHelpers.cpp
      main.cpp )

ADD_EXECUTABLE ( gl-nox.out
                 ${ginzburg_landau_SRCS} )

SET ( internal_LIBS
      gl
      cartesian
      io
    )

# ==============================================================================
FIND_PACKAGE( Boost COMPONENTS system filesystem REQUIRED )
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIRS} )
SET(Boost_USE_STATIC_LIBS   ON)
SET(Boost_USE_MULTITHREADED ON)
SET( LINK_LIBRARIES ${LINK_LIBRARIES}
                            ${Boost_LIBRARIES} )
# ------------------------------------------------------------------------------
FIND_LIBRARY( GFORTRAN ${GFORTRAN_LIBRARY_NAMES}
              HINTS ${GFORTRAN_LIBRARY_DIRS}    )
SET( LINK_LIBRARIES ${LINK_LIBRARIES}
                    ${GFORTRAN} )
# ==============================================================================
# Set executable linking information.
TARGET_LINK_LIBRARIES ( gl-nox.out
                        ${internal_LIBS}
                        LOCA
                        NOXEpetra
                        ${LINK_LIBRARIES}
                      )
# ==============================================================================
