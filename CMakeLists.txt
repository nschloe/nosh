# We need CONFIGURE_PACKAGE_CONFIG_FILE, available from CMake 2.8.8 on.
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)

# There's no Fortran in this project, but Trilinos has some Fortran in it
# (e.g., via BLAS/LAPACK).
# This is a cause for undefined symbols at the exectuable linking stage, mostly
# for Fortran's standard library (if the linker language is C++). Adding
# Fortran to the list of languages in this project makes sure that CMake
# correctly appends the needed libraries.
# This only influences static builds.
# At some point in time, the BLAS/LAPACK dependency on Fortran should be
# propagated up to applications like Nosh such that explicitly adding the
# Fortran compiler here isn't necessary anymore.
PROJECT(Nosh CXX)

IF (NOT Trilinos_Implicit)
  #FIND_PACKAGE(Trilinos REQUIRED)
  FIND_PACKAGE(Trilinos REQUIRED COMPONENTS Belos MueLu Thyra Tpetra NOX Piro)
ENDIF()
FIND_PACKAGE(Mikado REQUIRED)

FIND_PACKAGE(Eigen3 REQUIRED)
FIND_PACKAGE(MOAB REQUIRED)

# Make sure the compilers match.
IF(NOT ${Trilinos_CXX_COMPILER} STREQUAL ${CMAKE_CXX_COMPILER})
  MESSAGE(WARNING "C++ compilers don't match (Trilinos: ${Trilinos_CXX_COMPILER}, ${PROJECT_NAME}: ${CMAKE_CXX_COMPILER}).")
ENDIF()
#IF(NOT ${Trilinos_Fortran_COMPILER} STREQUAL ${CMAKE_Fortran_COMPILER})
#  MESSAGE(FATAL_ERROR "Fortran compilers don't match (Trilinos: ${Trilinos_Fortran_COMPILER}, ${PROJECT_NAME}: ${CMAKE_Fortran_COMPILER}).")
#ENDIF()

set(Nosh_MAJOR_VERSION 11)
set(Nosh_MINOR_VERSION 0)
set(Nosh_PATCH_VERSION 0)
set(
  Nosh_VERSION
  ${Nosh_MAJOR_VERSION}.${Nosh_MINOR_VERSION}.${Nosh_PATCH_VERSION}
  )

## Add compiler flags (e.g., -openmp).
#LIST(APPEND CMAKE_CXX_FLAGS ${Trilinos_CXX_COMPILER_FLAGS})

IF(NOT DEFINED CMAKE_BUILD_TYPE)
   SET(
     CMAKE_BUILD_TYPE Debug
     CACHE STRING "Choose the type of build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Profile Release RelWithDebInfo MinSizeRel.")
ENDIF()

ENABLE_TESTING()

IF(${TIME_MONITOR})
    # Add support for timing for non-debug builds.
    ADD_DEFINITIONS(-DNOSH_TEUCHOS_TIME_MONITOR)
ENDIF()

IF(CMAKE_COMPILER_IS_GNUCXX)
  #SET(CMAKE_CXX_FLAGS_DEBUG "-Og -g -ggdb -Wall -pedantic -fbounds-check -Wextra -Wstrict-null-sentinel -Wshadow -Woverloaded-virtual -Weffc++ -Wsign-compare -ansi -std=c++11" )
    SET(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -ggdb -Wall -pedantic -fbounds-check -Wextra -Wstrict-null-sentinel -Wshadow -Woverloaded-virtual -Weffc++ -Wsign-compare -ansi -std=c++11")
    execute_process(COMMAND ${CMAKE_C_COMPILER} -dumpversion
                    OUTPUT_VARIABLE GCC_VERSION
                    )
    if (GCC_VERSION VERSION_GREATER 4.8)
      SET(CMAKE_CXX_FLAGS_DEBUG ${CMAKE_CXX_FLAGS_DEBUG} -fdiagnostics-color=auto)
    endif()
    SET(CMAKE_CXX_FLAGS_PROFILE "-pg -O2" )
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
ELSEIF("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
    SET(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -ggdb -Wall -Wextra -pedantic")
ENDIF()

option(ENABLE_GCOV "Enable gcov support" OFF)
if(ENABLE_GCOV)
  message(STATUS "Enabling gcov support")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} --coverage")
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} --coverage")
endif()

# Offer the user the choice of overriding the installation directories
SET(INSTALL_LIB_DIR     "lib"          CACHE PATH "Installation directory for libraries")
SET(INSTALL_BIN_DIR     "bin"          CACHE PATH "Installation directory for executables")
SET(INSTALL_INCLUDE_DIR "include/nosh" CACHE PATH "Installation directory for header files")
if(WIN32 AND NOT CYGWIN)
  set(DEF_INSTALL_CMAKE_DIR CMake)
else()
  set(DEF_INSTALL_CMAKE_DIR "lib/cmake/${PROJECT_NAME}")
endif()
set(
  INSTALL_CMAKE_DIR ${DEF_INSTALL_CMAKE_DIR} CACHE PATH
  "Installation directory for CMake files"
  )

# Add subdirectories after INSTALL_BIN_DIR has been set.
ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(examples)
ADD_SUBDIRECTORY(executables)
ADD_SUBDIRECTORY(test)

# Create and install some helper files.
INCLUDE(CMakePackageConfigHelpers)
CONFIGURE_PACKAGE_CONFIG_FILE(
  "${PROJECT_NAME}Config.cmake.in"
  "${PROJECT_NAME}Config.cmake"
  INSTALL_DESTINATION "${PROJECT_BINARY_DIR}"
  PATH_VARS INSTALL_INCLUDE_DIR INSTALL_LIB_DIR
  )

WRITE_BASIC_PACKAGE_VERSION_FILE(
  ${CMAKE_CURRENT_BINARY_DIR}/NoshConfigVersion.cmake
  VERSION ${Nosh_VERSION}
  COMPATIBILITY SameMajorVersion
  )

INSTALL(
  FILES
  "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
  "${PROJECT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
  DESTINATION "${INSTALL_CMAKE_DIR}"
  )
