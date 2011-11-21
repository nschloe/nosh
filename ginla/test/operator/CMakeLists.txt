FIND_PACKAGE( Trilinos REQUIRED )
INCLUDE_DIRECTORIES ( ${Trilinos_INCLUDE_DIRS}
                      ${Trilinos_TPL_INCLUDE_DIRS} )
LINK_DIRECTORIES( ${Trilinos_LIBRARY_DIRS}
                  ${Trilinos_TPL_LIBRARY_DIRS} )

INCLUDE_DIRECTORIES ( ${ginzburg-landau_SOURCE_DIR}/src )

SET( OPERATORTEST_EXECUTABLE
     "operatorTest.exe" )

ADD_EXECUTABLE ( ${OPERATORTEST_EXECUTABLE}
                 Ginla_KeoUnitTests.cpp
                 Ginla_JacUnitTests.cpp
                 unitTestMain.cpp )

# copy over the binary data
if( NOT ${CMAKE_CURRENT_SOURCE_DIR} EQUAL ${CMAKE_CURRENT_BINARY_DIR} )
    FILE( COPY "${CMAKE_CURRENT_SOURCE_DIR}/data/cubesmall.e"
          DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
    FILE( COPY "${CMAKE_CURRENT_SOURCE_DIR}/data/cubesmall-balanced.par.2.0"
          DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
    FILE( COPY "${CMAKE_CURRENT_SOURCE_DIR}/data/cubesmall-balanced.par.2.1"
          DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
endif()

# add tests
ADD_TEST( operatorTest
          ${OPERATORTEST_EXECUTABLE}
        )
ADD_TEST( operatorTestMpi2
          mpiexec -n 2 ${OPERATORTEST_EXECUTABLE}
        )

SET ( internal_LIBS
      ginla
    )

# Set executable linking information.
TARGET_LINK_LIBRARIES ( ${OPERATORTEST_EXECUTABLE}
                        ${internal_LIBS}
                        ${Trilinos_LIBRARIES}
                        ${Trilinos_TPL_LIBRARIES}
                      )
