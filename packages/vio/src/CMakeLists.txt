SET ( vio_SRCS
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

#ADD_LIBRARY ( vio STATIC
ADD_LIBRARY ( vio
              ${vio_SRCS} )

# Get the absolute paths of the libraries.
# -- this is recommended for TARGET_LINK_LIBRARIES.
FIND_LIBRARY( VTKIO_LIBRARY vtkIO ${VTK_LIBRARY_DIRS} )
FIND_LIBRARY( VTKFILTERING_LIBRARY vtkFiltering ${VTK_LIBRARY_DIRS} )
FIND_LIBRARY( VTKCOMMON_LIBRARY vtkCommon ${VTK_LIBRARY_DIRS} )
FIND_LIBRARY( VTKIMAGING_LIBRARY vtkImaging ${VTK_LIBRARY_DIRS} )

# list the libraries against which libio must be linked here
TARGET_LINK_LIBRARIES ( vio
                        ${VTKIO_LIBRARY}
                        ${VTKFILTERING_LIBRARY}
                        ${VTKCOMMON_LIBRARY}
                        ${VTKIMAGING_LIBRARY} )
