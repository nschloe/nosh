CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

FIND_PACKAGE( Trilinos REQUIRED )

# set compilers *before the PROJECT is declared
SET( CMAKE_C_COMPILER ${Trilinos_C_COMPILER} )
SET( CMAKE_CXX_COMPILER ${Trilinos_CXX_COMPILER} )
SET( CMAKE_Fortran_COMPILER ${Trilinos_Fortran_COMPILER} )

PROJECT( ginzburg-landau CXX Fortran )

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

# Add preprocessor definition "_DEBUG_" for debug builds.
STRING(COMPARE EQUAL ${CMAKE_BUILD_TYPE} "Debug" IS_DEBUG_BUILD )
IF( IS_DEBUG_BUILD )
    ADD_DEFINITIONS( -D_DEBUG_ )
ENDIF()

IF(CMAKE_COMPILER_IS_GNUCXX)
    #SET(CMAKE_CXX_FLAGS_RELEASE     "-O2" )
    SET(CMAKE_CXX_FLAGS_DEBUG     "-g -ggdb -Wall -pedantic -fbounds-check -Wextra -Wno-unused -Wstrict-null-sentinel -Wshadow -Woverloaded-virtual -ansi -std=c++0x" )
    SET(CMAKE_CXX_FLAGS_PROFILE     "-pg -O2" )
    SET(CMAKE_CXX_FLAGS_DEBUG     "-g -ggdb -Wall -pedantic -Wextra -Wno-unused-parameter -Wno-unused -Wshadow -Woverloaded-virtual -ansi -std=c++0x" )
ENDIF()

#SET(CMAKE_Fortran_FLAGS_RELEASE "-O2" )
SET(CMAKE_Fortran_FLAGS_DEBUG "-g -ggdb -Wall -pedantic -fbounds-check -Wextra -Wno-unused" )
SET(CMAKE_Fortran_FLAGS_PROFILE "-pg -O2" )

add_subdirectory ( src )
add_subdirectory ( examples )
add_subdirectory ( test )
