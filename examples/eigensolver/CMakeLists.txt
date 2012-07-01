CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

INCLUDE_DIRECTORIES ( ${cuantico_SOURCE_DIR}/src
                      ${Trilinos_INCLUDE_DIRS}
                      ${Trilinos_TPL_INCLUDE_DIRS} )

LINK_DIRECTORIES( ${Trilinos_LIBRARY_DIRS}
                  ${Trilinos_TPL_LIBRARY_DIRS} )
# ------------------------------------------------------------------------------
SET( MY_EXECUTABLE "eigensolver" )
ADD_EXECUTABLE( ${MY_EXECUTABLE}
                "eigensolvers.cpp" )
TARGET_LINK_LIBRARIES( ${MY_EXECUTABLE}
                       "cuantico"
                       ${Trilinos_LIBRARIES}
                       ${Trilinos_TPL_LIBRARIES} )
# ------------------------------------------------------------------------------
