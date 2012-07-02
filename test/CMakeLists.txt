FIND_PACKAGE( Trilinos REQUIRED )
INCLUDE_DIRECTORIES ( ${cuantico_SOURCE_DIR}/src
                      ${Trilinos_INCLUDE_DIRS}
                      ${Trilinos_TPL_INCLUDE_DIRS} )
LINK_DIRECTORIES( ${Trilinos_LIBRARY_DIRS}
                  ${Trilinos_TPL_LIBRARY_DIRS} )

SET(internal_LIBS "cuantico")

# copy over the binary data
if( NOT ${CMAKE_CURRENT_SOURCE_DIR} EQUAL ${CMAKE_CURRENT_BINARY_DIR} )
    FILE( COPY "${CMAKE_CURRENT_SOURCE_DIR}/data/rectanglesmall.e"
          DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
    FILE( COPY "${CMAKE_CURRENT_SOURCE_DIR}/data/rectanglesmall-balanced.par.2.0"
          DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
    FILE( COPY "${CMAKE_CURRENT_SOURCE_DIR}/data/rectanglesmall-balanced.par.2.1"
          DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
    FILE( COPY "${CMAKE_CURRENT_SOURCE_DIR}/data/pacman.e"
          DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
    FILE( COPY "${CMAKE_CURRENT_SOURCE_DIR}/data/pacman-balanced.par.2.0"
          DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
    FILE( COPY "${CMAKE_CURRENT_SOURCE_DIR}/data/pacman-balanced.par.2.1"
          DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
    FILE( COPY "${CMAKE_CURRENT_SOURCE_DIR}/data/cubesmall.e"
          DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
    FILE( COPY "${CMAKE_CURRENT_SOURCE_DIR}/data/cubesmall-balanced.par.2.0"
          DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
    FILE( COPY "${CMAKE_CURRENT_SOURCE_DIR}/data/cubesmall-balanced.par.2.1"
          DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
    FILE( COPY "${CMAKE_CURRENT_SOURCE_DIR}/data/brick-w-hole.e"
          DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
    FILE( COPY "${CMAKE_CURRENT_SOURCE_DIR}/data/brick-w-hole-balanced.par.2.0"
          DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
    FILE( COPY "${CMAKE_CURRENT_SOURCE_DIR}/data/brick-w-hole-balanced.par.2.1"
          DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
endif()
# ------------------------------------------------------------------------------
SET( IOTEST_EXECUTABLE
     "ioTest" )
ADD_EXECUTABLE ( ${IOTEST_EXECUTABLE}
                 Cuantico_IOTests.cpp
                 unitTestMain.cpp )
# add tests
ADD_TEST( ioTestSerial
          ${IOTEST_EXECUTABLE}
        )
ADD_TEST( ioTestMpi2
          mpiexec -n 2 ${IOTEST_EXECUTABLE}
        )

# Set executable linking information.
TARGET_LINK_LIBRARIES ( ${IOTEST_EXECUTABLE}
                        ${internal_LIBS}
                        ${Trilinos_LIBRARIES}
                        ${Trilinos_TPL_LIBRARIES}
                      )
# ------------------------------------------------------------------------------
SET( MESHTEST_EXECUTABLE
     "meshTest" )
ADD_EXECUTABLE ( ${MESHTEST_EXECUTABLE}
                 Cuantico_StkMeshTests.cpp
#                  Cuantico_EdgeCacheTests.cpp
                 unitTestMain.cpp )
# add tests
ADD_TEST( meshTestSerial
          ${MESHTEST_EXECUTABLE}
        )
ADD_TEST( meshTestMpi2
          mpiexec -n 2 ${MESHTEST_EXECUTABLE}
        )

# Set executable linking information.
TARGET_LINK_LIBRARIES ( ${MESHTEST_EXECUTABLE}
                        ${internal_LIBS}
                        ${Trilinos_LIBRARIES}
                        ${Trilinos_TPL_LIBRARIES}
                      )
# ------------------------------------------------------------------------------
SET( OPERATORTEST_EXECUTABLE
     "operatorTest" )
ADD_EXECUTABLE ( ${OPERATORTEST_EXECUTABLE}
                 Cuantico_KeoUnitTests.cpp
                 Cuantico_ComputeFUnitTests.cpp
                 Cuantico_JacUnitTests.cpp
                 Cuantico_dfdpUnitTests.cpp
                 unitTestMain.cpp )

# add tests
ADD_TEST( operatorTest
          ${OPERATORTEST_EXECUTABLE}
        )
ADD_TEST( operatorTestMpi2
          mpiexec -n 2 ${OPERATORTEST_EXECUTABLE}
        )

# Set executable linking information.
TARGET_LINK_LIBRARIES ( ${OPERATORTEST_EXECUTABLE}
                        ${internal_LIBS}
                        ${Trilinos_LIBRARIES}
                        ${Trilinos_TPL_LIBRARIES}
                      )
# ------------------------------------------------------------------------------
