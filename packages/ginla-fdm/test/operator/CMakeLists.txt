FIND_PACKAGE( Trilinos REQUIRED )
INCLUDE_DIRECTORIES ( ${Trilinos_INCLUDE_DIRS}
                      ${Trilinos_TPL_INCLUDE_DIRS} )
LINK_DIRECTORIES( ${Trilinos_LIBRARY_DIRS}
                  ${Trilinos_TPL_LIBRARY_DIRS} )

FIND_PACKAGE( Boost
              COMPONENTS system filesystem unit_test_framework
              REQUIRED )
INCLUDE_DIRECTORIES ( ${Boost_INCLUDE_DIRS} )

INCLUDE_DIRECTORIES ( ${ginzburg-landau_SOURCE_DIR}/packages/ginla-fdm/src
                      ${ginzburg-landau_SOURCE_DIR}/packages/ginla-fdm/src-recti
                    )

SET( MY_EXECUTABLE
     "operatorTest.out" )

# see http://www.docstoc.com/docs/14857835/CTest-and-Dashboard-_Slides-by-J
#INCLUDE( CTest )
ADD_TEST( operatorTest
          ${MY_EXECUTABLE}
          "--xml-input-file=${CMAKE_CURRENT_SOURCE_DIR}/data/conf.xml"
        )

ADD_EXECUTABLE ( ${MY_EXECUTABLE}
                 main.cpp )

SET ( internal_LIBS
      ginlafdm
      komplex2
      recti
      vio
    )

# Set executable linking information.
TARGET_LINK_LIBRARIES ( ${MY_EXECUTABLE}
                        ${internal_LIBS}
                        ${Trilinos_LIBRARIES}
                        ${Trilinos_TPL_LIBRARIES}
                        ${Boost_LIBRARIES}
                        bfd
                      )
