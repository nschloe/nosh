CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

FIND_PACKAGE( NOX REQUIRED )
FIND_PACKAGE( VTK REQUIRED )

FIND_PACKAGE( Boost COMPONENTS system filesystem REQUIRED )
INCLUDE_DIRECTORIES ( ${Boost_INCLUDE_DIRS} )

INCLUDE_DIRECTORIES ( ${ginzburg-landau_SOURCE_DIR}/packages/vio/src
                      ${ginzburg-landau_SOURCE_DIR}/packages/komplex2/src
                      ${ginzburg-landau_SOURCE_DIR}/packages/recti/src
                      ${ginzburg-landau_SOURCE_DIR}/packages/ginla/src
                      ${ginzburg-landau_SOURCE_DIR}/packages/ginla/src-fdm
                      ${ginzburg-landau_SOURCE_DIR}/packages/ginla/src-fvm
                      ${VTK_INCLUDE_DIRS}
                      ${Trilinos_INCLUDE_DIRS}
                      ${Trilinos_TPL_INCLUDE_DIRS} )

LINK_DIRECTORIES( ${Boost_LIBRARY_DIRS}
                  ${VTK_LIBRARY_DIRS}
                  ${Trilinos_LIBRARY_DIRS}
                  ${Trilinos_TPL_LIBRARY_DIRS} )
# ------------------------------------------------------------------------------
SET ( pirodriver_SRCS
      piro-driver.cpp  )
SET ( pirodriver_LIBS
      ginlafdm
      ginla
      komplex2
      recti
      vio
    )
SET( MY_EXECUTABLE
     "piro-driver.exe" )
ADD_EXECUTABLE( ${MY_EXECUTABLE}
                ${pirodriver_SRCS} )
SET_TARGET_PROPERTIES( ${MY_EXECUTABLE}
                       PROPERTIES LINKER_LANGUAGE Fortran )
TARGET_LINK_LIBRARIES( ${MY_EXECUTABLE}
                       ${pirodriver_LIBS}
                       ${Trilinos_LIBRARIES}
                       ${Trilinos_TPL_LIBRARIES}
                       ${Boost_LIBRARIES} )
# ------------------------------------------------------------------------------
SET ( pirodriverfvm_SRCS
      piro-driver-fvm.cpp  )
SET ( pirodriverfvm_LIBS
      ginlafvm
      ginla
      komplex2
      recti
      vio
    )
SET( MY_EXECUTABLE
     "piro-driver-fvm.exe" )
ADD_EXECUTABLE( ${MY_EXECUTABLE}
                ${pirodriverfvm_SRCS} )
SET_TARGET_PROPERTIES( ${MY_EXECUTABLE}
                       PROPERTIES LINKER_LANGUAGE Fortran )
TARGET_LINK_LIBRARIES( ${MY_EXECUTABLE}
                       ${pirodriverfvm_LIBS}
                       ${Trilinos_LIBRARIES}
                       ${Trilinos_TPL_LIBRARIES}
                       ${Boost_LIBRARIES} )
# ------------------------------------------------------------------------------
