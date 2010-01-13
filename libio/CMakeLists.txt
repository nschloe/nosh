INCLUDE( ../cmake-system.txt REQUIRED )

ENABLE_LANGUAGE( CXX )

SET ( gl_SRCS
      ioVirtual.cpp
      ioVtk.cpp
      IoVtkStructuredPoints.cpp
      ioVti.cpp
      ioXdmf.cpp
      ioFactory.cpp )

FIND_PACKAGE( Trilinos REQUIRED )
INCLUDE_DIRECTORIES( ${Trilinos_INCLUDE_DIR} )

FIND_PACKAGE( Boost REQUIRED )
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIRS} )

FIND_PACKAGE( VTK REQUIRED )
INCLUDE_DIRECTORIES( ${VTK_INCLUDE_DIRS} )

FIND_PACKAGE( HDF5 REQUIRED )
#INCLUDE_DIRECTORIES( ${HDF5_DIR}/include )

IF(MPI_COMPILED)
    FIND_PACKAGE( MPI REQUIRED )
    INCLUDE_DIRECTORIES( ${MPI_INCLUDE_PATH} )
ENDIF(MPI_COMPILED)

ADD_LIBRARY ( io STATIC
              ${gl_SRCS} )

# list the libraries against which libio must be linked here
TARGET_LINK_LIBRARIES ( io 
                        hdf5 vtkIO vtkFiltering vtkCommon vtkImaging )
