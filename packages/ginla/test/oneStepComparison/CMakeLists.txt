FIND_PACKAGE( Trilinos REQUIRED )
INCLUDE_DIRECTORIES ( ${Trilinos_INCLUDE_DIRS}
                      ${Trilinos_TPL_INCLUDE_DIRS} )
LINK_DIRECTORIES( ${Trilinos_LIBRARY_DIRS}
                  ${Trilinos_TPL_LIBRARY_DIRS} )

FIND_PACKAGE( Boost COMPONENTS system filesystem REQUIRED )
SET(Boost_USE_STATIC_LIBS   ON)
SET(Boost_USE_MULTITHREADED ON)
INCLUDE_DIRECTORIES ( ${Boost_INCLUDE_DIRS} )

FIND_PACKAGE( VTK REQUIRED )
INCLUDE_DIRECTORIES( ${VTK_INCLUDE_DIRS} )

FIND_PACKAGE( Boost COMPONENTS unit_test_framework REQUIRED )

INCLUDE_DIRECTORIES ( ${ginzburg-landau_SOURCE_DIR}/packages/vio/src
                      ${ginzburg-landau_SOURCE_DIR}/packages/recti/src
                      ${ginzburg-landau_SOURCE_DIR}/packages/ginla/src )

SET( MY_EXECUTABLE
     "glTestSuite.out" )

# see http://www.docstoc.com/docs/14857835/CTest-and-Dashboard-_Slides-by-J
#INCLUDE( CTest )
ADD_TEST( glTestSuite
          ${MY_EXECUTABLE}
          "--xml-input-file=${CMAKE_CURRENT_SOURCE_DIR}/data/conf.xml"
          "--expected-solution-file=${CMAKE_CURRENT_SOURCE_DIR}/data/expected-solution.vti"
        )

ADD_EXECUTABLE ( ${MY_EXECUTABLE}
                 main.cpp )

SET ( internal_LIBS
      ginla
      recti
      vio
    )

SET_TARGET_PROPERTIES( ${MY_EXECUTABLE}
                       PROPERTIES LINKER_LANGUAGE Fortran )

# Set executable linking information.
TARGET_LINK_LIBRARIES ( ${MY_EXECUTABLE}
                        ${internal_LIBS}
                        ${Trilinos_LIBRARIES}
                        ${Trilinos_TPL_LIBRARIES}
                        ${Boost_LIBRARIES}
                      )
