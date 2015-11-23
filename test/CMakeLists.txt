FIND_PACKAGE(MPI REQUIRED)

INCLUDE_DIRECTORIES(${Nosh_SOURCE_DIR}/src/)
INCLUDE_DIRECTORIES(SYSTEM
  ${Trilinos_INCLUDE_DIRS}
  ${Trilinos_TPL_INCLUDE_DIRS}
  ${EIGEN3_INCLUDE_DIR}
  )

SET(internal_LIBS "nosh")

#set(input_files
#  rectanglesmall
#  pacman
#  screw
#  shell
#  sphere
#  cubesmall
#  brick-w-hole
#  )
#
## Download the files
#INCLUDE(ExternalData)
#set(
#  ExternalData_URL_TEMPLATES
#  "https://downloads.sourceforge.net/project/noshdata/%(algo)/%(hash)"
#  )
#foreach(file ${input_files})
#  ExternalData_Expand_Arguments(
#    noshTestFetchData
#    OUT_DATA DATA{${CMAKE_SOURCE_DIR}/test/data/${file}.e}
#    )
#endforeach()
#ExternalData_Add_Target(noshTestFetchData)
#
## convert and split the downloaded files
#add_custom_target(create_h5m ALL)
#add_dependencies(create_h5m noshTestFetchData)
#
#SET(FILES2)
#SET(FILES7)
#foreach(file ${input_files})
#  add_custom_command(
#    TARGET create_h5m
#    WORKING_DIRECTORY data/
#    COMMAND meshio-convert ${file}.e ${file}.h5m
#    )
#
#  LIST(APPEND FILES2 ${file}-2.h5m)
#  add_custom_command(
#    OUTPUT data/${file}-2.h5m
#    DEPENDS data/${file}.h5m
#    WORKING_DIRECTORY data/
#    COMMAND mbpart 2 -m ML_KWAY ${file}.h5m ${file}-2.h5m
#    )
#
#  LIST(APPEND FILES7 ${file}-7.h5m)
#  add_custom_command(
#    OUTPUT data/${file}-7.h5m
#    DEPENDS data/${file}.h5m
#    WORKING_DIRECTORY data/
#    COMMAND mbpart 7 -m ML_KWAY ${file}.h5m ${file}-7.h5m
#    )
#endforeach()
#
#add_custom_target(split2 ALL
#  # DEPENDS noshTestFetchData
#  SOURCES ${FILES2}
#  )
#
#add_custom_target(split7 ALL
#  #  DEPENDS noshTestFetchData
#  SOURCES ${FILES7}
#  )


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
SET(OPERATORTEST_EXECUTABLE "operatorTest")
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


ADD_SUBDIRECTORY(data)
