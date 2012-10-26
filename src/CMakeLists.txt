# Treat Trilinos as system code to prevent warnings from
# Trilinos to appear here.
INCLUDE_DIRECTORIES(SYSTEM
                    ${Trilinos_INCLUDE_DIRS})

FILE(GLOB nosh_SRCS "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")
FILE(GLOB nosh_HEADERS "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp")

ADD_LIBRARY(nosh
            ${nosh_SRCS})

# Add target link libraries.
# This has no effect for static builds, but appends the Trilinos
# libraries to the link line to avoid underlinking.
# Not adding TARGET_LINK_LIBRARIES() makes the link step fail with
# -DCMAKE_SHARED_LINKER_FLAGS:STRING="-Wl,--no-undefined".
# We don't need the TPL libraries here since the linking information
# about them is already included in the ${Trilinos_LIBRARIES}.
TARGET_LINK_LIBRARIES(nosh
                      ${Trilinos_LIBRARIES}
                     )

SET_TARGET_PROPERTIES(nosh PROPERTIES
  PUBLIC_HEADER "${nosh_HEADERS}"
  VERSION ${${PROJECT_NAME}_VERSION}
  SOVERSION ${${PROJECT_NAME}_MAJOR_VERSION})

# Installation of the program
INSTALL(TARGETS nosh
  # IMPORTANT: Add an the export (and install it below). It will contain
  # important build information.
  EXPORT NoshLibraryDepends
  RUNTIME DESTINATION "${INSTALL_BIN_DIR}" COMPONENT bin
  LIBRARY DESTINATION "${INSTALL_LIB_DIR}" COMPONENT shlib
  PUBLIC_HEADER DESTINATION "${INSTALL_INCLUDE_DIR}" COMPONENT dev
  #PUBLIC_HEADER DESTINATION "${INSTALL_INCLUDE_DIR}/nosh" COMPONENT dev
  )

# Install the export set.
INSTALL(EXPORT NoshLibraryDepends
        DESTINATION "${INSTALL_CMAKE_DIR}"
        )
