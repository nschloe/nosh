SET ( vio_SRCS
      VIO_Reader_Abstract.cpp
      VIO_Image_Writer_Factory.cpp
      VIO_Image_Writer_Abstract.cpp
      VIO_Image_Writer_VTI.cpp
      VIO_Image_Writer_legacyVTK.cpp
      VIO_Image_Reader_Factory.cpp
      VIO_Image_Reader_Abstract.cpp
      VIO_Image_Reader_VTI.cpp
      VIO_Image_Reader_legacyVTK.cpp
      VIO_Mesh_Mesh.cpp
      VIO_Mesh_Reader.cpp
      VIO_Mesh_Writer.cpp
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
FIND_LIBRARY( VTKGRAPHICS_LIBRARY vtkGraphics ${VTK_LIBRARY_DIRS} ) # for vtkFeatureEdges

# list the libraries against which libio must be linked here
TARGET_LINK_LIBRARIES ( vio
                        ${VTKIO_LIBRARY}
                        ${VTKFILTERING_LIBRARY}
                        ${VTKCOMMON_LIBRARY}
                        ${VTKIMAGING_LIBRARY}
                        ${VTKGRAPHICS_LIBRARY}
                      )
