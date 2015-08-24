FIND_PACKAGE(MPI REQUIRED)

INCLUDE_DIRECTORIES(${Nosh_SOURCE_DIR}/src/)
INCLUDE_DIRECTORIES(SYSTEM
  ${Trilinos_INCLUDE_DIRS}
  ${Trilinos_TPL_INCLUDE_DIRS}
  ${EIGEN3_INCLUDE_DIR}
  )

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

set(input_files
  rectanglesmall.e
  pacman.e
  shell.e
  sphere.e
  cubesmall.e
  brick-w-hole.e
  )

# Download the files
INCLUDE(ExternalData)
set(
  ExternalData_URL_TEMPLATES
  "https://downloads.sourceforge.net/project/noshdata/%(algo)/%(hash)"
  )
foreach(file ${input_files})
  ExternalData_Expand_Arguments(
    noshTestFetchData
    OUT_DATA DATA{${CMAKE_SOURCE_DIR}/test/data/${file}}
    )
endforeach()
ExternalData_Add_Target(noshTestFetchData)

# ------------------------------------------------------------------------------
SET(IOTEST_EXECUTABLE
  "ioTest")
ADD_EXECUTABLE (${IOTEST_EXECUTABLE}
  io.cpp
  unitTestMain.cpp
  )
# add tests
ADD_TEST(ioTestSerial
  ${IOTEST_EXECUTABLE}
  )
# Check out FindMPI.cmake for proper syntax.
ADD_TEST(ioTestMpi2
  ${Trilinos_MPI_EXEC} --noprefix ${Trilinos_MPI_EXEC_NUMPROCS_FLAG} 2 ${IOTEST_EXECUTABLE}
  )
ADD_TEST(ioTestMpi7
  ${Trilinos_MPI_EXEC} --noprefix ${Trilinos_MPI_EXEC_NUMPROCS_FLAG} 7 ${IOTEST_EXECUTABLE}
  )

# Set executable linking information.
TARGET_LINK_LIBRARIES(
  ${IOTEST_EXECUTABLE}
  ${internal_LIBS}
  )
IF (NOT Trilinos_Implicit)
  TARGET_LINK_LIBRARIES(
    ${IOTEST_EXECUTABLE}
    ${Trilinos_LIBRARIES}
    )
ENDIF()
# ------------------------------------------------------------------------------
SET(MESHTEST_EXECUTABLE
  "meshTest")
ADD_EXECUTABLE (${MESHTEST_EXECUTABLE}
  mesh.cpp
  #                 nosh_Edge_cacheTests.cpp
  unitTestMain.cpp)
# add tests
ADD_TEST(meshTestSerial
  ${MESHTEST_EXECUTABLE}
  )
ADD_TEST(meshTestMpi2
  ${Trilinos_MPI_EXEC} --noprefix ${Trilinos_MPI_EXEC_NUMPROCS_FLAG} 2 ${MESHTEST_EXECUTABLE}
  )
ADD_TEST(meshTestMpi7
  ${Trilinos_MPI_EXEC} --noprefix ${Trilinos_MPI_EXEC_NUMPROCS_FLAG} 7 ${MESHTEST_EXECUTABLE}
  )

# Set executable linking information.
TARGET_LINK_LIBRARIES(
  ${MESHTEST_EXECUTABLE}
  ${internal_LIBS}
  )
IF (NOT Trilinos_Implicit)
  TARGET_LINK_LIBRARIES(
    ${MESHTEST_EXECUTABLE}
    ${Trilinos_LIBRARIES}
    )
ENDIF()
# ------------------------------------------------------------------------------
SET(OPERATORTEST_EXECUTABLE
  "operatorTest")
ADD_EXECUTABLE(${OPERATORTEST_EXECUTABLE}
  keo.cpp
  compute_f.cpp
  jac.cpp
  dfdp.cpp
  unitTestMain.cpp
  )

# add tests
ADD_TEST(operatorTest
  ${OPERATORTEST_EXECUTABLE}
  )
ADD_TEST(operatorTestMpi2
  ${Trilinos_MPI_EXEC} --noprefix ${Trilinos_MPI_EXEC_NUMPROCS_FLAG} 2 ${OPERATORTEST_EXECUTABLE}
  )
ADD_TEST(operatorTestMpi7
  ${Trilinos_MPI_EXEC} --noprefix ${Trilinos_MPI_EXEC_NUMPROCS_FLAG} 7 ${OPERATORTEST_EXECUTABLE}
  )

# Set executable linking information.
TARGET_LINK_LIBRARIES(
  ${OPERATORTEST_EXECUTABLE}
  ${internal_LIBS}
  )
IF (NOT Trilinos_Implicit)
  TARGET_LINK_LIBRARIES(
    ${OPERATORTEST_EXECUTABLE}
    ${Trilinos_LIBRARIES}
    )
ENDIF()
# ------------------------------------------------------------------------------
