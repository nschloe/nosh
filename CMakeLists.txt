cmake_minimum_required( VERSION 2.6 )

# set compiler variables *before declaring the project
INCLUDE( SetTrilinosVariables )

PROJECT( ginzburg-landau )

#SET(CMAKE_BUILD_TYPE Profile)

SET(CMAKE_BUILD_TYPE Debug)
SET(CMAKE_CXX_FLAGS_DEBUG "-g -pg -ggdb -Wall" )
# SET(CMAKE_CXX_FLAGS_DEBUG "-g -pg -ggdb -Wall -fbounds-check -pedantic -Wextra -Wno-unused -std=c++0x -Weffc++ -Wstrict-null-sentinel" )
#SET(CMAKE_CXX_FLAGS_DEBUG "" )

add_subdirectory ( packages )
add_subdirectory ( executables )
