CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

INCLUDE_DIRECTORIES ( ${nosh_SOURCE_DIR}/src
                      ${Trilinos_INCLUDE_DIRS}
                      ${Trilinos_TPL_INCLUDE_DIRS} )

LINK_DIRECTORIES( ${Boost_LIBRARY_DIRS}
                  ${Trilinos_LIBRARY_DIRS}
                  ${Trilinos_TPL_LIBRARY_DIRS} )
# ------------------------------------------------------------------------------
SET(MY_EXECUTABLE "cont")
ADD_EXECUTABLE( ${MY_EXECUTABLE} "continuation.cpp" "MyScalarField.cpp")
TARGET_LINK_LIBRARIES( ${MY_EXECUTABLE}
                       "nosh"
                       ${Trilinos_LIBRARIES}
                       ${Trilinos_TPL_LIBRARIES} )
# ------------------------------------------------------------------------------
