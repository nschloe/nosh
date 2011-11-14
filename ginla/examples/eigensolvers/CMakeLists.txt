CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

INCLUDE_DIRECTORIES ( ${ginzburg-landau_SOURCE_DIR}/src
                      ${Trilinos_INCLUDE_DIRS}
                      ${Trilinos_TPL_INCLUDE_DIRS} )

LINK_DIRECTORIES( ${Trilinos_LIBRARY_DIRS}
                  ${Trilinos_TPL_LIBRARY_DIRS} )
# ------------------------------------------------------------------------------
SET ( eigensolver_SRCS
      eigensolvers.cpp  )
SET ( eigensolver_LIBS
      ginla
    )
SET( MY_EXECUTABLE
     "eigensolvers.exe" )
ADD_EXECUTABLE( ${MY_EXECUTABLE}
                ${eigensolver_SRCS} )
TARGET_LINK_LIBRARIES( ${MY_EXECUTABLE}
                       ${eigensolver_LIBS}
                       ${Trilinos_LIBRARIES}
                       ${Trilinos_TPL_LIBRARIES} )
# ------------------------------------------------------------------------------
