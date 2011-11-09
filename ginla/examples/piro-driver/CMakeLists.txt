CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

INCLUDE_DIRECTORIES ( ${ginzburg-landau_SOURCE_DIR}/src
                      ${Trilinos_INCLUDE_DIRS}
                      ${Trilinos_TPL_INCLUDE_DIRS} )

LINK_DIRECTORIES( ${Boost_LIBRARY_DIRS}
                  ${Trilinos_LIBRARY_DIRS}
                  ${Trilinos_TPL_LIBRARY_DIRS} )
# ------------------------------------------------------------------------------
SET ( pirodriver_SRCS
      piro-driver.cpp  )
SET ( pirodriver_LIBS
      ginla
    )
SET( MY_EXECUTABLE
     "piro-driver.exe" )
ADD_EXECUTABLE( ${MY_EXECUTABLE}
                ${pirodriver_SRCS} )
TARGET_LINK_LIBRARIES( ${MY_EXECUTABLE}
                       ${pirodriver_LIBS}
                       ${Trilinos_LIBRARIES}
                       ${Trilinos_TPL_LIBRARIES} )
# ------------------------------------------------------------------------------
