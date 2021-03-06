FIND_PACKAGE( Trilinos REQUIRED )
INCLUDE_DIRECTORIES ( ${Trilinos_INCLUDE_DIRS}
                      ${Trilinos_TPL_INCLUDE_DIRS} )
LINK_DIRECTORIES( ${Trilinos_LIBRARY_DIRS}
                  ${Trilinos_TPL_LIBRARY_DIRS} )

FIND_PACKAGE( Boost
              COMPONENTS system filesystem unit_test_framework
              REQUIRED )
INCLUDE_DIRECTORIES ( ${Boost_INCLUDE_DIRS} )

INCLUDE_DIRECTORIES ( ${ginzburg-landau_SOURCE_DIR}/packages/vio/src
                      ${ginzburg-landau_SOURCE_DIR}/packages/komplex2/src
                      ${ginzburg-landau_SOURCE_DIR}/packages/recti/src
                      ${ginzburg-landau_SOURCE_DIR}/packages/ginla/src
                      ${ginzburg-landau_SOURCE_DIR}/packages/ginla/src-fdm
                    )

SET( MY_EXECUTABLE
     "branchSwitch.out" )

# see http://www.docstoc.com/docs/14857835/CTest-and-Dashboard-_Slides-by-J
#INCLUDE( CTest )
ADD_TEST( branchSwitch
          ${MY_EXECUTABLE}
          "--xml-input-file=${CMAKE_CURRENT_SOURCE_DIR}/data/conf.xml"
          "--expected-solution-file=${CMAKE_CURRENT_SOURCE_DIR}/data/expected-solution.vti"
        )

ADD_EXECUTABLE ( ${MY_EXECUTABLE}
                 main.cpp )

SET ( internal_LIBS
      recti
      ginlafdm
      ginla
      komplex2
      vio
    )

# Set executable linking information.
TARGET_LINK_LIBRARIES ( ${MY_EXECUTABLE}
                        ${internal_LIBS}
                        ${Trilinos_LIBRARIES}
                        ${Trilinos_TPL_LIBRARIES}
                        ${Boost_LIBRARIES}
                      )
