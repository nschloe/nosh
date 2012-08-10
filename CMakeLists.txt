CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

FIND_PACKAGE(Trilinos REQUIRED)

# Set compilers before the PROJECT is declared.
SET(CMAKE_CXX_COMPILER ${Trilinos_CXX_COMPILER})

PROJECT(Nosh CXX)

set(NOSH_MAJOR_VERSION 10)
set(NOSH_MINOR_VERSION 0)
set(NOSH_PATCH_VERSION 0)
set(Nosh_VERSION
  ${NOSH_MAJOR_VERSION}.${NOSH_MINOR_VERSION}.${NOSH_PATCH_VERSION})

ENABLE_TESTING()

IF(DEFINED CMAKE_BUILD_TYPE)
   SET(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING "Choose the type of
build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug
Profile Release RelWithDebInfo MinSizeRel.")
ELSE()
   SET(CMAKE_BUILD_TYPE Debug CACHE STRING "Choose the type of build,
options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Profile
Release RelWithDebInfo MinSizeRel.")
ENDIF()

STRING(COMPARE EQUAL ${CMAKE_BUILD_TYPE} "Debug" IS_DEBUG_BUILD)
IF( IS_DEBUG_BUILD )
    # Add preprocessor definition "_DEBUG_" for debug builds.
    ADD_DEFINITIONS( -D_DEBUG_ )
    # Don't do optimization for debug builds.
    SET(CMAKE_CXX_FLAGS_DEBUG "-O0 ${CMAKE_CXX_FLAGS_DEBUG}")
ELSE()
    # Add support for timing for non-debug builds.
    ADD_DEFINITIONS( -DNOSH_TEUCHOS_TIME_MONITOR )
ENDIF()

IF(CMAKE_COMPILER_IS_GNUCXX)
    SET(CMAKE_CXX_FLAGS_DEBUG     "-O0 -g -ggdb -Wall -pedantic -fbounds-check -Wextra -Wno-unused -Wstrict-null-sentinel -Wshadow -Woverloaded-virtual -ansi -std=c++0x" )
    SET(CMAKE_CXX_FLAGS_PROFILE     "-pg -O2" )
ENDIF()

# Offer the user the choice of overriding the installation directories
set(INSTALL_LIB_DIR lib CACHE PATH "Installation directory for libraries")
set(INSTALL_BIN_DIR bin CACHE PATH "Installation directory for executables")
set(INSTALL_INCLUDE_DIR include CACHE PATH "Installation directory for header files")
if(WIN32 AND NOT CYGWIN)
  set(DEF_INSTALL_CMAKE_DIR CMake)
else()
  set(DEF_INSTALL_CMAKE_DIR "lib/cmake/Nosh")
endif()
set(INSTALL_CMAKE_DIR ${DEF_INSTALL_CMAKE_DIR} CACHE PATH
  "Installation directory for CMake files")

# Make relative paths absolute (needed later on)
foreach(p LIB BIN INCLUDE CMAKE)
  set(var INSTALL_${p}_DIR)
  if(NOT IS_ABSOLUTE "${${var}}")
    set(${var} "${CMAKE_INSTALL_PREFIX}/${${var}}")
  endif()
endforeach()

add_subdirectory( src )
add_subdirectory( examples )
add_subdirectory( test )


# Add all targets to the build-tree export set
export(TARGETS nosh cont
       FILE "${PROJECT_BINARY_DIR}/NoshLibraryDepends.cmake")
 
# Export the package for use from the build-tree
# (this registers the build-tree with a global CMake-registry)
export(PACKAGE nosh)
 
# Create a NoshBuildTreeSettings.cmake file for the use from the build tree
configure_file(NoshBuildTreeSettings.cmake.in
               "${PROJECT_BINARY_DIR}/NoshTreeSettings.cmake" @ONLY)

# Set some variables that are later exported
SET(Nosh_CONFIG_LIBRARIES "nosh")
SET(Nosh_CONFIG_Trilinos_INCLUDE_DIRS ${Trilinos_INCLUDE_DIRS})
SET(Nosh_CONFIG_Trilinos_LIBRARY_DIRS ${Trilinos_LIBRARY_DIRS})
SET(Nosh_CONFIG_Trilinos_LIBRARIES ${Trilinos_LIBRARIES})

# Create the NoshConfig.cmake and NoshConfigVersion files
file(RELATIVE_PATH Nosh_CONFIG_INCLUDE_DIRS
     "${INSTALL_CMAKE_DIR}"
     "${INSTALL_INCLUDE_DIR}")
file(RELATIVE_PATH Nosh_CONFIG_LIBRARY_DIRS
     "${INSTALL_CMAKE_DIR}"
     "${INSTALL_LIB_DIR}")
configure_file(NoshConfig.cmake.in
               "${PROJECT_BINARY_DIR}/NoshConfig.cmake")
configure_file(NoshConfigVersion.cmake.in
               "${PROJECT_BINARY_DIR}/NoshConfigVersion.cmake" @ONLY)

# Install the export set for use with the install-tree
install(EXPORT NoshLibraryDepends
        DESTINATION "lib/cmake/Nosh"
        COMPONENT dev)
        #"${INSTALL_DATA_DIR}/Nosh/CMake"

install(FILES
        "${PROJECT_BINARY_DIR}/NoshConfig.cmake"
        "${PROJECT_BINARY_DIR}/NoshConfigVersion.cmake"
        DESTINATION "${INSTALL_CMAKE_DIR}")
