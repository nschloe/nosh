CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

INCLUDE_DIRECTORIES ( ${Trilinos_INCLUDE_DIRS}
                      ${Trilinos_TPL_INCLUDE_DIRS} )

LINK_DIRECTORIES( ${Trilinos_LIBRARY_DIRS}
                  ${Trilinos_TPL_LIBRARY_DIRS} )
# ------------------------------------------------------------------------------
SET ( scaltest_SRCS
      dummy.f
      scalability-tests.cpp  )
SET( MY_EXECUTABLE
     "scaltest.exe" )
ADD_EXECUTABLE( ${MY_EXECUTABLE}
                ${scaltest_SRCS} )
TARGET_LINK_LIBRARIES( ${MY_EXECUTABLE}
                       ${Trilinos_LIBRARIES}
                       ${Trilinos_TPL_LIBRARIES}
                    )
# ------------------------------------------------------------------------------
