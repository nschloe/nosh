CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

FIND_PACKAGE( Trilinos REQUIRED )

# set compilers *before the PROJECT is declared
SET( CMAKE_C_COMPILER ${Trilinos_C_COMPILER} )
SET( CMAKE_CXX_COMPILER ${Trilinos_CXX_COMPILER} )
SET( CMAKE_Fortran_COMPILER ${Trilinos_Fortran_COMPILER} )

PROJECT( ginzburg-landau CXX Fortran )

ENABLE_TESTING()

SET(CMAKE_BUILD_TYPE Release)
# SET(CMAKE_BUILD_TYPE Debug)
# SET(CMAKE_BUILD_TYPE Profile)

SET(CMAKE_CXX_FLAGS_RELEASE     "-O2" )
SET(CMAKE_Fortran_FLAGS_RELEASE "-O2" )

SET(CMAKE_CXX_FLAGS_DEBUG     "-g -ggdb -Wall -pedantic -ansi -fbounds-check -Wextra -Wno-unused -Wstrict-null-sentinel -std=c++0x" )
SET(CMAKE_Fortran_FLAGS_DEBUG "-g -ggdb -Wall -pedantic -ansi -fbounds-check -Wextra -Wno-unused -Wstrict-null-sentinel" )

SET(CMAKE_CXX_FLAGS_PROFILE     "-pg -O2" )
SET(CMAKE_Fortran_FLAGS_PROFILE "-pg -O2" )

add_subdirectory ( example )
add_subdirectory ( src )
add_subdirectory ( src-komplex2 )
add_subdirectory ( src-recti )
add_subdirectory ( src-vio )
add_subdirectory ( test )
