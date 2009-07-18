ENABLE_LANGUAGE( CXX )

include_directories ( ${ginzburg-landau_SOURCE_DIR}/libgl )

SET ( ginzburg_landau_SRCS
      main.cpp )

add_executable ( gl.out
                 ${ginzburg_landau_SRCS} )

SET ( internal_LIBS
      gl
    )

target_link_libraries ( gl.out
			${internal_LIBS}
                      )
