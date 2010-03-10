FIND_PACKAGE( NOX REQUIRED )
INCLUDE_DIRECTORIES ( ${NOX_INCLUDE_DIRS}
                      ${NOX_TPL_INCLUDE_DIRS} )
LINK_DIRECTORIES( ${NOX_LIBRARY_DIRS}
                  ${NOX_TPL_LIBRARY_DIRS} )

FIND_PACKAGE( VTK REQUIRED )
INCLUDE_DIRECTORIES( ${VTK_INCLUDE_DIRS} )

FIND_PACKAGE( Boost COMPONENTS unit_test_framework REQUIRED )



# INCLUDE_DIRECTORIES ( ${ginzburg-landau_SOURCE_DIR}/packages/libio/src
#                       ${ginzburg-landau_SOURCE_DIR}/packages/libcartesian/src
#                       ${ginzburg-landau_SOURCE_DIR}/packages/libgl/src )

SET( MY_EXECUTABLE
     "glTestSuite.out" )

# see http://www.docstoc.com/docs/14857835/CTest-and-Dashboard-_Slides-by-J
#INCLUDE( CTest )
ADD_TEST( glTestSuite
          ${MY_EXECUTABLE} )

ADD_EXECUTABLE ( ${MY_EXECUTABLE}
                 main.cpp )

SET ( internal_LIBS )

SET_TARGET_PROPERTIES( ${MY_EXECUTABLE}
                       PROPERTIES LINKER_LANGUAGE Fortran )

# Set executable linking information.
TARGET_LINK_LIBRARIES ( ${MY_EXECUTABLE}
                        ${internal_LIBS}
                        ${NOX_LIBRARIES}
                        ${NOX_TPL_LIBRARIES}
                        ${Boost_LIBRARIES}
                      )
