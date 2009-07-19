ENABLE_LANGUAGE( CXX )

SET ( gl_SRCS 
      psiGrid.cpp 
      aGrid.cpp 
      ginzburgLandau.cpp
      glSystem.cpp )

INCLUDE_DIRECTORIES( /usr/local/trilinos/090603a/include/ )

add_library ( gl STATIC
              ${gl_SRCS} )

target_link_libraries ( gl )
