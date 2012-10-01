CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)

# There's no Fortran in this project, but Trilinos has some Fortran
# in it (e.g., via BLAS/LAPACK).
# This is a cause for undefined symbols at the exectuable linking
# stage, mostly for Fortran's standard library (if the linker language
# is C++). Adding Fortran to the list of languages in this project makes
# sure that CMake correctly appends the needed libraries.
# This only influences static builds.
# At some point in time, the BLAS/LAPACK dependency on Fortran should
# be propagated up to applications like Nosh such that explicitly
# adding the Fortran compiler here isn't necessary anymore.
PROJECT(Nosh CXX Fortran)

FIND_PACKAGE(Trilinos REQUIRED)

# Make sure the compilers match.
IF(NOT ${Trilinos_CXX_COMPILER} STREQUAL ${CMAKE_CXX_COMPILER})
  MESSAGE(FATAL_ERROR "C++ compilers don't match (Trilinos: ${Trilinos_CXX_COMPILER}, Nosh: ${CMAKE_CXX_COMPILER}).")
ENDIF()
IF(NOT ${Trilinos_Fortran_COMPILER} STREQUAL ${CMAKE_Fortran_COMPILER})
  MESSAGE(FATAL_ERROR "Fortran compilers don't match (Trilinos: ${Trilinos_Fortran_COMPILER}, Nosh: ${CMAKE_Fortran_COMPILER}).")
ENDIF()

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

ENABLE_TESTING()

#IF(NOT ${CMAKE_BUILD_TYPE} STREQUAL "Debug")
#    # Add support for timing for non-debug builds.
#    ADD_DEFINITIONS(-DNOSH_TEUCHOS_TIME_MONITOR)
#ENDIF()

IF(CMAKE_COMPILER_IS_GNUCXX)
    SET(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -ggdb -Wall -pedantic -fbounds-check -Wextra -Wstrict-null-sentinel -Wshadow -Woverloaded-virtual -Weffc++ -ansi -std=c++0x" )
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
