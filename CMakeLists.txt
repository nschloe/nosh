cmake_minimum_required( VERSION 2.6 )

PROJECT( ginzburg-landau )

# set compilers
#SET(CMAKE_CXX_COMPILER /opt/openmpi-1.3.3/bin/mpic++)

#SET(CMAKE_BUILD_TYPE Profile)

SET(CMAKE_BUILD_TYPE Debug)
#SET(CMAKE_CXX_FLAGS_DEBUG "-pg -ggdb -Wall -Wextra -Weffc++ -pedantic -fbounds-check" )
SET(CMAKE_CXX_FLAGS_DEBUG "-pg -ggdb -Wall -Wextra" )

add_subdirectory ( libgl )
add_subdirectory ( main )
add_subdirectory ( nox )
add_subdirectory ( loca )
