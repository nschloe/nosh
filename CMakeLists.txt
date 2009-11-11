cmake_minimum_required( VERSION 2.6 )

PROJECT( ginzburg-landau )

INCLUDE( cmake-system.txt REQUIRED )

IF(MPI_COMPILED)
    FIND_PACKAGE( MPI REQUIRED )
    #     SET(CMAKE_CXX_COMPILER ${MPI_COMPILER} )
    #    SET(CMAKE_CXX_COMPILER /usr/local/bin/mpic++)
ENDIF(MPI_COMPILED)


# set compilers
SET(CMAKE_CXX_COMPILER /opt/gcc/4.4.2/bin/g++ )
SET(CMAKE_C_COMPILER /opt/gcc/4.4.2/bin/gcc )

#SET(CMAKE_BUILD_TYPE Profile)

SET(CMAKE_BUILD_TYPE Debug)
#SET(CMAKE_CXX_FLAGS_DEBUG "-pg -ggdb -Wall -Wextra -Weffc++ -pedantic -fbounds-check" )
SET(CMAKE_CXX_FLAGS_DEBUG "-g -pg -ggdb -Wall -Wextra" )
#SET(CMAKE_CXX_FLAGS_DEBUG "" )

SET(GL_CMAKE_DIR "${ginzburg-landau_SOURCE_DIR}/CMakeModules")
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
