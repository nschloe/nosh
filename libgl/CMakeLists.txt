ENABLE_LANGUAGE( CXX )

SET ( gl_SRCS 
      psiGrid.cpp 
      aGrid.cpp 
      ginzburgLandau.cpp
      glSystem.cpp )

INCLUDE_DIRECTORIES( /usr/local/trilinos/090603a/include/ )

ADD_LIBRARY ( gl STATIC
              ${gl_SRCS} )

TARGET_LINK_LIBRARIES ( gl )
