CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

FIND_PACKAGE( Boost COMPONENTS system filesystem REQUIRED )
INCLUDE_DIRECTORIES ( ${Boost_INCLUDE_DIRS} )

INCLUDE_DIRECTORIES ( ${Trilinos_INCLUDE_DIRS}
                      ${Trilinos_TPL_INCLUDE_DIRS} )

LINK_DIRECTORIES( ${Boost_LIBRARY_DIRS}
                  ${VTK_LIBRARY_DIRS}
                  ${Trilinos_LIBRARY_DIRS}
                  ${Trilinos_TPL_LIBRARY_DIRS} )
# ------------------------------------------------------------------------------
SET( MY_EXECUTABLE
     "cg-indefinite-test.exe" )
ADD_EXECUTABLE( ${MY_EXECUTABLE}
                cg-indefinite-test.cpp )
TARGET_LINK_LIBRARIES( ${MY_EXECUTABLE}
                       ${Trilinos_LIBRARIES}
                       ${Trilinos_TPL_LIBRARIES}
                       ${Boost_LIBRARIES} )
# ------------------------------------------------------------------------------
