INCLUDE( ${ginzburg-landau_SOURCE_DIR}/cmake-system.txt REQUIRED )

ENABLE_LANGUAGE( CXX )

SET ( cartesian_SRCS
      GridVirtual.cpp
      GridUniformVirtual.cpp
      GridSquare.cpp
      GridUniformSquare.cpp
      GridReader.cpp )

INCLUDE( SetTrilinosVariables REQUIRED )
INCLUDE_DIRECTORIES( ${Trilinos_INCLUDE_DIR} )

INCLUDE_DIRECTORIES( ${ginzburg-landau_SOURCE_DIR}/packages/libio/src )

ADD_LIBRARY ( cartesian STATIC
              ${cartesian_SRCS} )
