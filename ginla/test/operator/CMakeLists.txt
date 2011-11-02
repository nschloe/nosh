FIND_PACKAGE( Trilinos REQUIRED )
INCLUDE_DIRECTORIES ( ${Trilinos_INCLUDE_DIRS}
                      ${Trilinos_TPL_INCLUDE_DIRS} )
LINK_DIRECTORIES( ${Trilinos_LIBRARY_DIRS}
                  ${Trilinos_TPL_LIBRARY_DIRS} )

FIND_PACKAGE( Boost
              COMPONENTS unit_test_framework
              REQUIRED )
INCLUDE_DIRECTORIES ( ${Boost_INCLUDE_DIRS} )

INCLUDE_DIRECTORIES ( ${ginzburg-landau_SOURCE_DIR}/src )

SET( MY_EXECUTABLE
     "keoTest.exe" )

ADD_EXECUTABLE ( ${MY_EXECUTABLE}
                 main.cpp )

# see http://www.docstoc.com/docs/14857835/CTest-and-Dashboard-_Slides-by-J
ADD_TEST( keoTest
          ${MY_EXECUTABLE}
          --input=${CMAKE_CURRENT_SOURCE_DIR}/data/cubesmall.e
        )

ADD_TEST( keoTestMpi2
          mpiexec -n 2 ${MY_EXECUTABLE}
          --input=${CMAKE_CURRENT_SOURCE_DIR}/data/cubesmall-balanced.par
        )

SET ( internal_LIBS
      ginla
    )

# Set executable linking information.
TARGET_LINK_LIBRARIES ( ${MY_EXECUTABLE}
                        ${internal_LIBS}
                        ${Trilinos_LIBRARIES}
                        ${Trilinos_TPL_LIBRARIES}
                        ${Boost_LIBRARIES}
                      )
