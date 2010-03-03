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

FIND_PACKAGE( Trilinos REQUIRED )
INCLUDE_DIRECTORIES( ${Trilinos_INCLUDE_DIRS} )

FIND_PACKAGE( Boost REQUIRED )
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIRS} )

FIND_PACKAGE( VTK REQUIRED )
INCLUDE_DIRECTORIES( ${VTK_INCLUDE_DIRS} )
# TODO Use absolute paths in TARGET_LINK_LIBRARIES rather than LINK_LIBRARIES

ADD_LIBRARY ( io STATIC
              ${gl_SRCS} )

# list the libraries against which libio must be linked here
TARGET_LINK_LIBRARIES ( io 
                        vtkIO vtkFiltering vtkCommon vtkImaging )
