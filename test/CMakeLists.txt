FIND_PACKAGE(MPI REQUIRED)

INCLUDE_DIRECTORIES(${Nosh_SOURCE_DIR}/src/)
INCLUDE_DIRECTORIES(SYSTEM
  ${Trilinos_INCLUDE_DIRS}
  ${Trilinos_TPL_INCLUDE_DIRS}
  ${EIGEN3_INCLUDE_DIR}
  )

SET(internal_LIBS "nosh")


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

SET(FTEST_EXECUTABLE "fTest")
ADD_EXECUTABLE(${FTEST_EXECUTABLE}
  compute_f.cpp
  unitTestMain.cpp
  )
# Set executable linking information.
TARGET_LINK_LIBRARIES(
  ${FTEST_EXECUTABLE}
  ${internal_LIBS}
  )
IF (NOT Trilinos_Implicit)
  TARGET_LINK_LIBRARIES(
    ${FTEST_EXECUTABLE}
    ${Trilinos_LIBRARIES}
    )
ENDIF()
# add tests
ADD_TEST(fTest
  ${FTEST_EXECUTABLE}
  )
ADD_TEST(fTestMpi2
  ${Trilinos_MPI_EXEC} --noprefix ${Trilinos_MPI_EXEC_NUMPROCS_FLAG} 2 ${FTEST_EXECUTABLE}
  )
ADD_TEST(fTestMpi7
  ${Trilinos_MPI_EXEC} --noprefix ${Trilinos_MPI_EXEC_NUMPROCS_FLAG} 7 ${FTEST_EXECUTABLE}
  )


SET(KEOTEST_EXECUTABLE "keoTest")
ADD_EXECUTABLE(${KEOTEST_EXECUTABLE}
  keo.cpp
  unitTestMain.cpp
  )
# Set executable linking information.
TARGET_LINK_LIBRARIES(
  ${KEOTEST_EXECUTABLE}
  ${internal_LIBS}
  )
IF (NOT Trilinos_Implicit)
  TARGET_LINK_LIBRARIES(
    ${KEOTEST_EXECUTABLE}
    ${Trilinos_LIBRARIES}
    )
ENDIF()
# add tests
ADD_TEST(operatorTest
  ${KEOTEST_EXECUTABLE}
  )
ADD_TEST(operatorTestMpi2
  ${Trilinos_MPI_EXEC} --noprefix ${Trilinos_MPI_EXEC_NUMPROCS_FLAG} 2 ${KEOTEST_EXECUTABLE}
  )
ADD_TEST(operatorTestMpi7
  ${Trilinos_MPI_EXEC} --noprefix ${Trilinos_MPI_EXEC_NUMPROCS_FLAG} 7 ${KEOTEST_EXECUTABLE}
  )

SET(OPERATORTEST_EXECUTABLE "operatorTest")
ADD_EXECUTABLE(${OPERATORTEST_EXECUTABLE}
  jac.cpp
  dfdp.cpp
  unitTestMain.cpp
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


ADD_SUBDIRECTORY(data)
