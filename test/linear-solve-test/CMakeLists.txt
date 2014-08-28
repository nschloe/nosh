CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

INCLUDE_DIRECTORIES ( ${nosh_SOURCE_DIR}/src
  ${Trilinos_INCLUDE_DIRS}
  ${Trilinos_TPL_INCLUDE_DIRS} )

LINK_DIRECTORIES( ${Trilinos_LIBRARY_DIRS}
  ${Trilinos_TPL_LIBRARY_DIRS} )
# ------------------------------------------------------------------------------
SET( MY_EXECUTABLE
  "keo-aztecoo.exe" )
ADD_EXECUTABLE( ${MY_EXECUTABLE}
  "keo-aztecoo.cpp" )
TARGET_LINK_LIBRARIES( ${MY_EXECUTABLE}
  "nosh"
  ${Trilinos_LIBRARIES}
  ${Trilinos_TPL_LIBRARIES} )
# ------------------------------------------------------------------------------
SET( MY_EXECUTABLE
  "linop-tester.exe" )
ADD_EXECUTABLE( ${MY_EXECUTABLE}
  "linop-tester.cpp" )
TARGET_LINK_LIBRARIES( ${MY_EXECUTABLE}
  "nosh"
  ${Trilinos_LIBRARIES}
  ${Trilinos_TPL_LIBRARIES} )
# ------------------------------------------------------------------------------
SET( MY_EXECUTABLE
  "poisson1d-tester.exe" )
ADD_EXECUTABLE( ${MY_EXECUTABLE}
  "poisson1d-tester.cpp" )
TARGET_LINK_LIBRARIES( ${MY_EXECUTABLE}
  ${Trilinos_LIBRARIES}
  ${Trilinos_TPL_LIBRARIES} )
# ------------------------------------------------------------------------------
