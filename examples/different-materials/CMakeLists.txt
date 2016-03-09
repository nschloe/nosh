INCLUDE_DIRECTORIES(
  ${PROJECT_SOURCE_DIR}/src/
  ${CMAKE_CURRENT_BINARY_DIR}  # for the output header
  )
INCLUDE_DIRECTORIES(
  SYSTEM
  ${Trilinos_INCLUDE_DIRS}
  ${Trilinos_TPL_INCLUDE_DIRS}
  ${EIGEN3_INCLUDE_DIR}
  )

# generate header file
ADD_CUSTOM_COMMAND(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/different-materials.hpp
  COMMAND ${PROJECT_SOURCE_DIR}/nfc/nfc ${CMAKE_CURRENT_SOURCE_DIR}/different-materials.py ${CMAKE_CURRENT_BINARY_DIR}/different-materials.hpp
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/different-materials.py
  )

# copy over mesh files
ADD_CUSTOM_COMMAND(
   OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/pacman.e
   COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_BINARY_DIR}/test/data/pacman.e ${CMAKE_CURRENT_BINARY_DIR}/pacman.e
   DEPENDS noshTestFetchData
   )

SET(MY_EXECUTABLE different-materials)
ADD_EXECUTABLE(${MY_EXECUTABLE} different-materials.cpp different-materials.hpp pacman.e)
TARGET_LINK_LIBRARIES(${MY_EXECUTABLE} nosh)
