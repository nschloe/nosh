ENABLE_LANGUAGE( CXX )

SET ( gl_SRCS
      WriterFactory.cpp
      AbstractImageWriter.cpp
      VtiWriter.cpp
      VtkWriter.cpp
      ReaderFactory.cpp
      AbstractImageReader.cpp
      VtiReader.cpp
      VtkReader.cpp
    )

#       ioVirtual.cpp
#       ioXdmf.cpp
#       ioFactory.cpp

INCLUDE( SetTrilinosVariables REQUIRED )
INCLUDE_DIRECTORIES( ${Trilinos_INCLUDE_DIR} )

FIND_PACKAGE( Boost REQUIRED )
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIRS} )

FIND_PACKAGE( VTK REQUIRED )
INCLUDE_DIRECTORIES( ${VTK_INCLUDE_DIRS} )

ADD_LIBRARY ( io STATIC
              ${gl_SRCS} )

# list the libraries against which libio must be linked here
TARGET_LINK_LIBRARIES ( io 
                        vtkIO vtkFiltering vtkCommon vtkImaging )
