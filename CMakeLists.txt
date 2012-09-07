CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)

# Explicitly set this property to TRUE.
# This should be the case by default, but a CMake bug
# doesn't set it. (2012-09-04)
# If this is fixed in a future version of CMake,
# this following statement can be removed.
SET_PROPERTY(GLOBAL PROPERTY
             TARGET_SUPPORTS_SHARED_LIBS TRUE)

FIND_PACKAGE(Trilinos REQUIRED)

# Set compilers before the PROJECT is declared.
SET(CMAKE_CXX_COMPILER ${Trilinos_CXX_COMPILER})

set(Nosh_MAJOR_VERSION 11)
set(Nosh_MINOR_VERSION 0)
set(Nosh_PATCH_VERSION 0)
set(Nosh_VERSION
  ${Nosh_MAJOR_VERSION}.${Nosh_MINOR_VERSION}.${Nosh_PATCH_VERSION})


IF(NOT DEFINED CMAKE_BUILD_TYPE)
   SET(CMAKE_BUILD_TYPE None CACHE STRING "Choose the type of build,
options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Profile
Release RelWithDebInfo MinSizeRel.")
ENDIF()

PROJECT(Nosh CXX)
ENABLE_TESTING()

#IF(NOT ${CMAKE_BUILD_TYPE} STREQUAL "Debug")
#    # Add support for timing for non-debug builds.
#    ADD_DEFINITIONS(-DNOSH_TEUCHOS_TIME_MONITOR)
#ENDIF()

IF(CMAKE_COMPILER_IS_GNUCXX)
    SET(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -ggdb -Wall -pedantic -fbounds-check -Wextra -Wstrict-null-sentinel -Wshadow -Woverloaded-virtual -Weffc++ -ansi -std=c++11" )
    SET(CMAKE_CXX_FLAGS_PROFILE "-pg -O2" )
ENDIF()

# Offer the user the choice of overriding the installation directories
SET(INSTALL_LIB_DIR lib CACHE PATH "Installation directory for libraries")
SET(INSTALL_BIN_DIR bin CACHE PATH "Installation directory for executables")
SET(INSTALL_INCLUDE_DIR include CACHE PATH "Installation directory for header files")
if(WIN32 AND NOT CYGWIN)
  set(DEF_INSTALL_CMAKE_DIR CMake)
else()
  set(DEF_INSTALL_CMAKE_DIR "lib/cmake/Nosh")
endif()
set(INSTALL_CMAKE_DIR ${DEF_INSTALL_CMAKE_DIR} CACHE PATH
  "Installation directory for CMake files")

# Add subdirectories after INSTALL_BIN_DIR has been set.
ADD_SUBDIRECTORY( src )
ADD_SUBDIRECTORY( executables )
ADD_SUBDIRECTORY( test )

# Add all targets to the build-tree export set
export(TARGETS nosh cont
       FILE "${PROJECT_BINARY_DIR}/NoshLibraryDepends.cmake")
 
# Export the package for use from the build-tree
# (this registers the build-tree with a global CMake-registry)
export(PACKAGE nosh)

# Set some variables that are later exported
SET(Nosh_CONFIG_LIBRARIES "nosh")

INCLUDE(CMakePackageConfigHelpers)
CONFIGURE_PACKAGE_CONFIG_FILE("NoshConfig.cmake.in"
                              "NoshConfig.cmake"
                              INSTALL_DESTINATION "${PROJECT_BINARY_DIR}"
                              PATH_VARS INSTALL_INCLUDE_DIR INSTALL_LIB_DIR)

WRITE_BASIC_PACKAGE_VERSION_FILE(${CMAKE_CURRENT_BINARY_DIR}/NoshConfigVersion.cmake
                                 VERSION ${Nosh_VERSION}
                                 COMPATIBILITY SameMajorVersion)

# Install the export set for use with the install-tree.
INSTALL(EXPORT NoshLibraryDepends
        DESTINATION "lib/cmake/Nosh"
        COMPONENT dev)

INSTALL(FILES
        "${PROJECT_BINARY_DIR}/NoshConfig.cmake"
        "${PROJECT_BINARY_DIR}/NoshConfigVersion.cmake"
        DESTINATION "${INSTALL_CMAKE_DIR}")
