INCLUDE( ../cmake-system.txt REQUIRED )

ENABLE_LANGUAGE( CXX )

SET ( gl_SRCS
      glException.cpp
      staggeredGrid.cpp
      ioVirtual.cpp
      ioVtk.cpp
      ioVti.cpp
      ioXdmf.cpp
      ioFactory.cpp
      ginzburgLandau.cpp
      glBoundaryConditionsVirtual.cpp
      glBoundaryConditionsInner.cpp
      glBoundaryConditionsOuter.cpp
      glBoundaryConditionsCentral.cpp
      glPrePostOperator.cpp
      glSystem.cpp
      domainRectangle.cpp
      eigenSaver.cpp          )

FIND_PACKAGE( Trilinos REQUIRED )
INCLUDE_DIRECTORIES( ${Trilinos_INCLUDE_DIR} )

FIND_PACKAGE( Boost REQUIRED )
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIRS} )

IF(MPI_COMPILED)
    FIND_PACKAGE( MPI REQUIRED )
    INCLUDE_DIRECTORIES( ${MPI_INCLUDE_PATH} )
ENDIF(MPI_COMPILED)

INCLUDE_DIRECTORIES( ${PARACONT_DIR}/include )
INCLUDE_DIRECTORIES( ${HDF5_DIR}/include )

ADD_LIBRARY ( gl STATIC
              ${gl_SRCS} )

TARGET_LINK_LIBRARIES ( gl )
