INCLUDE( ../cmake-system.txt REQUIRED )

ENABLE_LANGUAGE( CXX )

IF(MPI_COMPILED)
    FIND_PACKAGE( MPI REQUIRED )
    SET( LINK_LIBRARIES ${LINK_LIBRARIES}
                        ${MPI_LIBRARY}    )
    INCLUDE_DIRECTORIES ( ${MPI_INCLUDE_PATH} )
ENDIF(MPI_COMPILED)

FIND_PACKAGE( LOCA      REQUIRED )
FIND_PACKAGE( NOXEpetra REQUIRED )


INCLUDE_DIRECTORIES ( ${Trilinos_INCLUDE_DIR} )

INCLUDE_DIRECTORIES ( ${ginzburg-landau_SOURCE_DIR}/libio
                      ${ginzburg-landau_SOURCE_DIR}/libgl )

SET ( ginzburg_landau_SRCS
      glNox.cpp
      main.cpp )

ADD_EXECUTABLE ( gl-nox.out
                 ${ginzburg_landau_SRCS} )

SET ( internal_LIBS
      gl
      io
    )

# ==============================================================================
FIND_LIBRARY( GFORTRAN ${GFORTRAN_LIBNAME}
              HINTS ${GFORTRAN_DIR}    )
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
