CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

FIND_PACKAGE( Trilinos REQUIRED )

# set compilers *before the PROJECT is declared
SET( CMAKE_C_COMPILER ${Trilinos_C_COMPILER} )
SET( CMAKE_CXX_COMPILER ${Trilinos_CXX_COMPILER} )
SET( CMAKE_Fortran_COMPILER ${Trilinos_Fortran_COMPILER} )

PROJECT( ginzburg-landau CXX Fortran )

ENABLE_TESTING()

#SET(CMAKE_BUILD_TYPE Profile)

SET(CMAKE_BUILD_TYPE Debug)
SET(CMAKE_CXX_FLAGS_DEBUG "-g -pg -ggdb -Wall" )
# SET(CMAKE_CXX_FLAGS_DEBUG "-g -pg -ggdb -Wall -fbounds-check -pedantic -Wextra -Wno-unused -std=c++0x -Weffc++ -Wstrict-null-sentinel" )
#SET(CMAKE_CXX_FLAGS_DEBUG "" )

add_subdirectory ( packages )
add_subdirectory ( executables )
