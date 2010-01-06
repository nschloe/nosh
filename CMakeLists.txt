cmake_minimum_required( VERSION 2.6 )

PROJECT( ginzburg-landau )

INCLUDE( cmake-system.txt REQUIRED )

#SET(CMAKE_BUILD_TYPE Profile)

SET(CMAKE_BUILD_TYPE Debug)
#SET(CMAKE_CXX_FLAGS_DEBUG "-pg -ggdb -Wall -Wextra -Weffc++ -pedantic -fbounds-check" )
SET(CMAKE_CXX_FLAGS_DEBUG "-g -pg -ggdb -Wall -Wextra -Wno-unused" )
#SET(CMAKE_CXX_FLAGS_DEBUG "" )

SET(GL_CMAKE_DIR "${ginzburg-landau_SOURCE_DIR}/../CMakeMacros")
SET(CMAKE_MODULE_PATH "${GL_CMAKE_DIR}" ${CMAKE_MODULE_PATH})

add_subdirectory ( libio )
add_subdirectory ( libgl )
add_subdirectory ( nox )
#add_subdirectory ( pseudo-timestepping )
add_subdirectory ( oneParameterContinuation )
add_subdirectory ( turning-point-continuation )
#add_subdirectory ( pitchfork-continuation )
#add_subdirectory ( branch-switching )
#add_subdirectory ( paracont )
