SET ( cartesian_SRCS
      DomainVirtual.cpp
      DomainEllipse.cpp
      DomainCircle.cpp
      DomainPolygon.cpp
      DomainRectangle.cpp
      DomainSquare.cpp
      GridVirtual.cpp
      Grid.cpp
      GridUniform.cpp
      GridUniformVirtual.cpp
      GridSquare.cpp
      GridUniformSquare.cpp
      GridReader.cpp )

FIND_PACKAGE( Trilinos REQUIRED )
INCLUDE_DIRECTORIES( ${Trilinos_INCLUDE_DIRS} )

FIND_PACKAGE( VTK REQUIRED )
INCLUDE_DIRECTORIES( ${VTK_INCLUDE_DIRS} )

INCLUDE_DIRECTORIES( ${ginzburg-landau_SOURCE_DIR}/packages/libio/src )

ADD_LIBRARY ( cartesian STATIC
              ${cartesian_SRCS} )
