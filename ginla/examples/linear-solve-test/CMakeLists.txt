CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

INCLUDE_DIRECTORIES ( ${ginzburg-landau_SOURCE_DIR}/src
                      ${Trilinos_INCLUDE_DIRS}
                      ${Trilinos_TPL_INCLUDE_DIRS} )

LINK_DIRECTORIES( ${Trilinos_LIBRARY_DIRS}
                  ${Trilinos_TPL_LIBRARY_DIRS} )
# ------------------------------------------------------------------------------
SET ( keoaztecoo_SRCS
      keo-aztecoo.cpp  )
SET ( keoaztecoo_LIBS
      ginla
    )
SET( MY_EXECUTABLE
     "keo-aztecoo.exe" )
ADD_EXECUTABLE( ${MY_EXECUTABLE}
                ${keoaztecoo_SRCS} )
TARGET_LINK_LIBRARIES( ${MY_EXECUTABLE}
                       ${keoaztecoo_LIBS}
                       ${Trilinos_LIBRARIES}
                       ${Trilinos_TPL_LIBRARIES} )
# ------------------------------------------------------------------------------
SET ( linoptester_SRCS
      linop-tester.cpp  )
SET ( linoptester_LIBS
      ginla
    )
SET( MY_EXECUTABLE
     "linop-tester.exe" )
ADD_EXECUTABLE( ${MY_EXECUTABLE}
                ${linoptester_SRCS} )
TARGET_LINK_LIBRARIES( ${MY_EXECUTABLE}
                       ${linoptester_LIBS}
                       ${Trilinos_LIBRARIES}
                       ${Trilinos_TPL_LIBRARIES} )
# ------------------------------------------------------------------------------
