INCLUDE( ../cmake-system.txt REQUIRED )

ENABLE_LANGUAGE( CXX )

SET ( gl_SRCS
      ioVirtual.cpp
      ioVtk.cpp
      IoVtkStructuredPoints.cpp
      ioVti.cpp
      ioXdmf.cpp
      ioFactory.cpp )

INCLUDE( SetTrilinosVariables REQUIRED )
INCLUDE_DIRECTORIES( ${Trilinos_INCLUDE_DIR} )

FIND_PACKAGE( Boost REQUIRED )
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIRS} )

FIND_PACKAGE( VTK REQUIRED )
INCLUDE_DIRECTORIES( ${VTK_INCLUDE_DIRS} )

FIND_PACKAGE( HDF5 REQUIRED )
INCLUDE_DIRECTORIES( ${HDF5_INCLUDE_DIR} )

ADD_LIBRARY ( io STATIC
              ${gl_SRCS} )

# list the libraries against which libio must be linked here
TARGET_LINK_LIBRARIES ( io 
                        ${HDF5_CXX_LIBRARIES}
                        vtkIO vtkFiltering vtkCommon vtkImaging )
