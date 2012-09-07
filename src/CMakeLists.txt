FIND_PACKAGE(Trilinos REQUIRED)

# Treat Trilinos as system code to prevent warnings from
# Trilinos to appear here.
INCLUDE_DIRECTORIES(SYSTEM
                    ${Trilinos_INCLUDE_DIRS}
                    ${Trilinos_TPL_INCLUDE_DIRS})

FILE(GLOB nosh_SRCS "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")
FILE(GLOB nosh_HEADERS "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp")

ADD_LIBRARY(nosh
            ${nosh_SRCS})

SET_TARGET_PROPERTIES(
  nosh
  PROPERTIES
  PUBLIC_HEADER "${nosh_HEADERS}"
  VERSION ${${PROJECT_NAME}_VERSION}
  SOVERSION ${${PROJECT_NAME}_MAJOR_VERSION})

# Installation of the program
INSTALL(TARGETS nosh
  DESTINATION "lib"
  # IMPORTANT: Add the nosh library to the "export-set"
  EXPORT NoshLibraryDepends
  RUNTIME DESTINATION "${INSTALL_BIN_DIR}" COMPONENT bin
  LIBRARY DESTINATION "${INSTALL_LIB_DIR}" COMPONENT shlib
  PUBLIC_HEADER DESTINATION "${INSTALL_INCLUDE_DIR}" COMPONENT dev)
  #PUBLIC_HEADER DESTINATION "${INSTALL_INCLUDE_DIR}/nosh" COMPONENT dev)
