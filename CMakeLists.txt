CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

FIND_PACKAGE( Trilinos REQUIRED )

# Set compilers before the PROJECT is declared.
SET(CMAKE_CXX_COMPILER ${Trilinos_CXX_COMPILER})

PROJECT(nosh CXX)

ENABLE_TESTING()

IF(DEFINED CMAKE_BUILD_TYPE)
   SET(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING "Choose the type of
build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug
Release RelWithDebInfo MinSizeRel.")
ELSE()
   SET(CMAKE_BUILD_TYPE Debug CACHE STRING "Choose the type of build,
options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release
RelWithDebInfo MinSizeRel.")
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

add_subdirectory( src )
add_subdirectory( examples )
add_subdirectory( test )
