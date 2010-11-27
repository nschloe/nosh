SET ( vio_SRCS
      VIO_TpetraMesh_Mesh.cpp
      VIO_TpetraMesh_Reader.cpp
      VIO_TpetraMesh_Writer.cpp
    )

FIND_PACKAGE( Trilinos REQUIRED )
INCLUDE_DIRECTORIES( ${Trilinos_INCLUDE_DIRS} )

FIND_PACKAGE( Boost REQUIRED )
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIRS} )

FIND_PACKAGE( VTK REQUIRED )
INCLUDE_DIRECTORIES( ${VTK_INCLUDE_DIRS} )

#ADD_LIBRARY ( vio2 STATIC
ADD_LIBRARY ( vio2
              ${vio_SRCS} )

# Get the absolute paths of the libraries.
# -- this is recommended for TARGET_LINK_LIBRARIES.
FIND_LIBRARY( VTKIO_LIBRARY vtkIO ${VTK_LIBRARY_DIRS} )
FIND_LIBRARY( VTKFILTERING_LIBRARY vtkFiltering ${VTK_LIBRARY_DIRS} )
FIND_LIBRARY( VTKCOMMON_LIBRARY vtkCommon ${VTK_LIBRARY_DIRS} )
FIND_LIBRARY( VTKIMAGING_LIBRARY vtkImaging ${VTK_LIBRARY_DIRS} )
FIND_LIBRARY( VTKGRAPHICS_LIBRARY vtkGraphics ${VTK_LIBRARY_DIRS} ) # for vtkFeatureEdges

# list the libraries against which libio must be linked here
TARGET_LINK_LIBRARIES ( vio2
                        ${VTKIO_LIBRARY}
                        ${VTKFILTERING_LIBRARY}
                        ${VTKCOMMON_LIBRARY}
                        ${VTKIMAGING_LIBRARY}
                        ${VTKGRAPHICS_LIBRARY}
                      )
