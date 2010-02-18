cmake_minimum_required( VERSION 2.6 )

PROJECT( ginzburg-landau )

#SET(CMAKE_BUILD_TYPE Profile)

SET(CMAKE_BUILD_TYPE Debug)
SET(CMAKE_CXX_FLAGS_DEBUG "-g -pg -ggdb -Wall" )
# SET(CMAKE_CXX_FLAGS_DEBUG "-g -pg -ggdb -Wall -fbounds-check -pedantic -Wextra -Wno-unused -std=c++0x -Weffc++ -Wstrict-null-sentinel" )
#SET(CMAKE_CXX_FLAGS_DEBUG "" )

SET(GL_CMAKE_DIR "${ginzburg-landau_SOURCE_DIR}/../CMakeMacros")
SET(CMAKE_MODULE_PATH "${GL_CMAKE_DIR}" ${CMAKE_MODULE_PATH})

add_subdirectory ( packages/libio )
add_subdirectory ( packages/libcartesian )
add_subdirectory ( packages/libgl )
add_subdirectory ( executables/nox )
#add_subdirectory ( pseudo-timestepping )
add_subdirectory ( executables/oneParameterContinuation )
# add_subdirectory ( executables/nox-predictor-step )
add_subdirectory ( executables/turning-point-continuation )
#add_subdirectory ( pitchfork-continuation )
#add_subdirectory ( branch-switching )
