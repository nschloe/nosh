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

SET( KEOTEST_EXECUTABLE
     "keoTest.exe" )

SET( JACTEST_EXECUTABLE
     "jacTest.exe" )

ADD_EXECUTABLE ( ${KEOTEST_EXECUTABLE}
                 keoTest.cpp )

ADD_EXECUTABLE ( ${JACTEST_EXECUTABLE}
                 jacTest.cpp )

# see http://www.docstoc.com/docs/14857835/CTest-and-Dashboard-_Slides-by-J
ADD_TEST( keoTest
          ${KEOTEST_EXECUTABLE}
          --input=${CMAKE_CURRENT_SOURCE_DIR}/data/cubesmall.e
        )

ADD_TEST( keoTestMpi2
          mpiexec -n 2 ${KEOTEST_EXECUTABLE}
          --input=${CMAKE_CURRENT_SOURCE_DIR}/data/cubesmall-balanced.par
        )

ADD_TEST( jacTest
          ${JACTEST_EXECUTABLE}
          --input=${CMAKE_CURRENT_SOURCE_DIR}/data/cubesmall.e
        )

ADD_TEST( jacTestMpi2
          mpiexec -n 2 ${JACTEST_EXECUTABLE}
          --input=${CMAKE_CURRENT_SOURCE_DIR}/data/cubesmall-balanced.par
        )

SET ( internal_LIBS
      ginla
    )

# Set executable linking information.
TARGET_LINK_LIBRARIES ( ${KEOTEST_EXECUTABLE}
                        ${internal_LIBS}
                        ${Trilinos_LIBRARIES}
                        ${Trilinos_TPL_LIBRARIES}
                        ${Boost_LIBRARIES}
                      )

TARGET_LINK_LIBRARIES ( ${JACTEST_EXECUTABLE}
                        ${internal_LIBS}
                        ${Trilinos_LIBRARIES}
                        ${Trilinos_TPL_LIBRARIES}
                        ${Boost_LIBRARIES}
                      )
