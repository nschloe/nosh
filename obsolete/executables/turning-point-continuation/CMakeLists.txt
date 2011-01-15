CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

FIND_PACKAGE( NOX REQUIRED )
FIND_PACKAGE( VTK REQUIRED )

FIND_PACKAGE( Boost COMPONENTS system filesystem REQUIRED )
INCLUDE_DIRECTORIES ( ${Boost_INCLUDE_DIRS} )

SET ( internal_LIBS
      ginlafdm
      ginla
      recti
      vio
    )

INCLUDE_DIRECTORIES ( ${ginzburg-landau_SOURCE_DIR}/packages/vio/src
                      ${ginzburg-landau_SOURCE_DIR}/packages/recti/src
                      ${ginzburg-landau_SOURCE_DIR}/packages/ginla/src
                      ${ginzburg-landau_SOURCE_DIR}/packages/ginla/src-fdm
                      ${VTK_INCLUDE_DIRS}
                      ${Trilinos_INCLUDE_DIRS}
                      ${Trilinos_TPL_INCLUDE_DIRS} )

LINK_DIRECTORIES( ${Boost_LIBRARY_DIRS}
                  ${VTK_LIBRARY_DIRS}
                  ${Trilinos_LIBRARY_DIRS}
                  ${Trilinos_TPL_LIBRARY_DIRS} )

SET ( opc_SRCS
      glTurningPointContinuation.cpp  )
SET( MY_EXECUTABLE
     "glTurningPointContinuation.out" )
ADD_EXECUTABLE( ${MY_EXECUTABLE}
                ${opc_SRCS} )
TARGET_LINK_LIBRARIES( ${MY_EXECUTABLE}
                       ${internal_LIBS}
                       ${Trilinos_LIBRARIES}
                       ${Trilinos_TPL_LIBRARIES}
                       ${Boost_LIBRARIES} )

SET ( piro_SRCS
      glTurningPointContinuationPiro.cpp )
SET( MY_EXECUTABLE "glTurningPointContinuationPiro.out" )
ADD_EXECUTABLE( ${MY_EXECUTABLE}
                ${piro_SRCS})
TARGET_LINK_LIBRARIES( ${MY_EXECUTABLE}
                       ${internal_LIBS}
                       ${Trilinos_LIBRARIES}
                       ${Trilinos_TPL_LIBRARIES}
                       ${Boost_LIBRARIES} )
