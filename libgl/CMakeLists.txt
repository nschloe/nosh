ENABLE_LANGUAGE( CXX )

SET ( gl_SRCS 
      staggeredGrid.cpp
      stateFileReader.cpp
      stateFileWriter.cpp
      ginzburgLandau.cpp
      glSystem.cpp )

INCLUDE_DIRECTORIES( /usr/local/trilinos/090603a/include/
                     /opt/trilinos-090603a/LINUX_MPI_64BIT/include/
                     /opt/boost/1.39.0/include/boost-1_39/ )

ADD_LIBRARY ( gl STATIC
              ${gl_SRCS} )

TARGET_LINK_LIBRARIES ( gl )
