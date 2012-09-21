
INCLUDE_DIRECTORIES(${Nosh_SOURCE_DIR}/src)
INCLUDE_DIRECTORIES(SYSTEM
                    ${Trilinos_INCLUDE_DIRS}
                    ${Trilinos_TPL_INCLUDE_DIRS})

## Get the absolute path of Trilinos_LIBRARY_DIRS. This is to avoid
## conflicts with Trilinos_TPL_LIBRARY_DIRS which is typically given
## as absolute paths.
## To be removed once Sandia bug #5692, cf.
##   https://software.sandia.gov/bugzilla/show_bug.cgi?id=5692,
## has been resolved.
#GET_FILENAME_COMPONENT(Trilinos_LIBRARY_DIRS_ABSOLUTE
#                       ${Trilinos_LIBRARY_DIRS}
#                       REALPATH)
#
#LINK_DIRECTORIES(${Trilinos_LIBRARY_DIRS_ABSOLUTE}
#                 ${Trilinos_TPL_LIBRARY_DIRS})

SET(internal_LIBS "nosh")

# copy over the binary data
if( NOT ${CMAKE_CURRENT_SOURCE_DIR} EQUAL ${CMAKE_CURRENT_BINARY_DIR} )
    FILE( COPY "${CMAKE_CURRENT_SOURCE_DIR}/data/rectanglesmall.e"
          DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
    FILE( COPY "${CMAKE_CURRENT_SOURCE_DIR}/data/pacman.e"
          DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
    FILE( COPY "${CMAKE_CURRENT_SOURCE_DIR}/data/shell.e"
          DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
    FILE( COPY "${CMAKE_CURRENT_SOURCE_DIR}/data/sphere.e"
          DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
    FILE( COPY "${CMAKE_CURRENT_SOURCE_DIR}/data/cubesmall.e"
          DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
    FILE( COPY "${CMAKE_CURRENT_SOURCE_DIR}/data/brick-w-hole.e"
          DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
endif()
# ------------------------------------------------------------------------------
SET( IOTEST_EXECUTABLE
     "ioTest" )
ADD_EXECUTABLE ( ${IOTEST_EXECUTABLE}
                 Nosh_IOTests.cpp
                 unitTestMain.cpp )
# add tests
ADD_TEST( ioTestSerial
          ${IOTEST_EXECUTABLE}
        )
ADD_TEST( ioTestMpi2
          mpiexec -n 2 ${IOTEST_EXECUTABLE}
        )
ADD_TEST( ioTestMpi7
          mpiexec -n 7 ${IOTEST_EXECUTABLE}
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
                 Nosh_StkMeshTests.cpp
#                 Nosh_EdgeCacheTests.cpp
                 unitTestMain.cpp )
# add tests
ADD_TEST( meshTestSerial
          ${MESHTEST_EXECUTABLE}
        )
ADD_TEST( meshTestMpi2
          mpiexec -n 2 ${MESHTEST_EXECUTABLE}
        )
ADD_TEST( meshTestMpi7
          mpiexec -n 7 ${MESHTEST_EXECUTABLE}
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
                 Nosh_KeoUnitTests.cpp
                 Nosh_ComputeFUnitTests.cpp
                 Nosh_JacUnitTests.cpp
                 Nosh_dfdpUnitTests.cpp
                 unitTestMain.cpp )

# add tests
ADD_TEST( operatorTest
          ${OPERATORTEST_EXECUTABLE}
        )
ADD_TEST( operatorTestMpi2
          mpiexec -n 2 ${OPERATORTEST_EXECUTABLE}
        )
ADD_TEST( operatorTestMpi7
          mpiexec -n 7 ${OPERATORTEST_EXECUTABLE}
        )

# Set executable linking information.
TARGET_LINK_LIBRARIES ( ${OPERATORTEST_EXECUTABLE}
                        ${internal_LIBS}
                        ${Trilinos_LIBRARIES}
                        ${Trilinos_TPL_LIBRARIES}
                      )
# ------------------------------------------------------------------------------
