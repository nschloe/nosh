INCLUDE( ../cmake-system.txt REQUIRED )

ENABLE_LANGUAGE( CXX )

SET ( gl_SRCS
      ioException.cpp
      ioVirtual.cpp
      ioVtk.cpp
      ioVti.cpp
      ioXdmf.cpp
      ioFactory.cpp )

FIND_PACKAGE( Trilinos REQUIRED )
INCLUDE_DIRECTORIES( ${Trilinos_INCLUDE_DIR} )

FIND_PACKAGE( Boost REQUIRED )
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIRS} )

IF(MPI_COMPILED)
    FIND_PACKAGE( MPI REQUIRED )
    INCLUDE_DIRECTORIES( ${MPI_INCLUDE_PATH} )
ENDIF(MPI_COMPILED)

INCLUDE_DIRECTORIES( ${HDF5_DIR}/include )

ADD_LIBRARY ( io STATIC
              ${gl_SRCS} )

TARGET_LINK_LIBRARIES ( io )
