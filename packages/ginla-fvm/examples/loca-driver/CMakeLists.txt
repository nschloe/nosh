CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

FIND_PACKAGE( Boost COMPONENTS system filesystem REQUIRED )
INCLUDE_DIRECTORIES ( ${Boost_INCLUDE_DIRS} )

INCLUDE_DIRECTORIES ( ${ginzburg-landau_SOURCE_DIR}/packages/ginla-fvm/src
                      ${Trilinos_INCLUDE_DIRS}
                      ${Trilinos_TPL_INCLUDE_DIRS} )

LINK_DIRECTORIES( ${Boost_LIBRARY_DIRS}
                  ${Trilinos_LIBRARY_DIRS}
                  ${Trilinos_TPL_LIBRARY_DIRS} )
# ------------------------------------------------------------------------------
SET ( locadriver_SRCS
      loca-driver-fvm.cpp  )
SET ( locadriver_LIBS
      ginlafvm
    )
SET( MY_EXECUTABLE
     "loca-driver-fvm.exe" )
ADD_EXECUTABLE( ${MY_EXECUTABLE}
                ${locadriver_SRCS} )
#SET_TARGET_PROPERTIES( ${MY_EXECUTABLE}
#                       PROPERTIES LINKER_LANGUAGE Fortran )
TARGET_LINK_LIBRARIES( ${MY_EXECUTABLE}
                       ${locadriver_LIBS}
                       ${Trilinos_LIBRARIES}
                       ${Trilinos_TPL_LIBRARIES}
                       ${Boost_LIBRARIES} )
# ------------------------------------------------------------------------------
