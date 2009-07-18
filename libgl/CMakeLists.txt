ENABLE_LANGUAGE( CXX )

SET ( gl_SRCS 
      psiGrid.cpp 
      aGrid.cpp 
      ginzburgLandau.cpp )

add_library ( gl STATIC
              ${gl_SRCS} )

target_link_libraries ( gl )
