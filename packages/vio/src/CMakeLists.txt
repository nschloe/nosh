SET ( gl_SRCS
      VIO_Writer_Factory.cpp
      VIO_Writer_Abstract.cpp
      VIO_Writer_VTI.cpp
      VIO_Writer_legacyVTK.cpp
      VIO_Reader_Factory.cpp
      VIO_Reader_Abstract.cpp
      VIO_Reader_VTI.cpp
      VIO_Reader_legacyVTK.cpp
    )

FIND_PACKAGE( Trilinos REQUIRED )
INCLUDE_DIRECTORIES( ${Trilinos_INCLUDE_DIRS} )

FIND_PACKAGE( Boost REQUIRED )
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIRS} )

FIND_PACKAGE( VTK REQUIRED )
INCLUDE_DIRECTORIES( ${VTK_INCLUDE_DIRS} )
# TODO Use absolute paths in TARGET_LINK_LIBRARIES rather than LINK_LIBRARIES

ADD_LIBRARY ( vio STATIC
              ${gl_SRCS} )

# list the libraries against which libio must be linked here
TARGET_LINK_LIBRARIES ( vio
                        vtkIO vtkFiltering vtkCommon vtkImaging )
