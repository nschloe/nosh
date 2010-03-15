CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

FIND_PACKAGE( NOX REQUIRED )
LINK_DIRECTORIES( ${Trilinos_LIBRARY_DIRS}
                  ${Trilinos_TPL_LIBRARY_DIRS} )
INCLUDE_DIRECTORIES ( ${Trilinos_INCLUDE_DIRS}
                      ${Trilinos_TPL_INCLUDE_DIRS} )
                      
FIND_PACKAGE( Boost COMPONENTS system filesystem REQUIRED )
SET(Boost_USE_STATIC_LIBS   ON)
SET(Boost_USE_MULTITHREADED ON)
INCLUDE_DIRECTORIES ( ${Boost_INCLUDE_DIRS} )
LINK_DIRECTORIES( ${Boost_LIBRARY_DIRS} )


SET ( internal_LIBS
      ginla
      recti
      vio
    )

INCLUDE_DIRECTORIES ( ${ginzburg-landau_SOURCE_DIR}/packages/vio/src
                      ${ginzburg-landau_SOURCE_DIR}/packages/recti/src
                      ${ginzburg-landau_SOURCE_DIR}/packages/ginla/src )

SET ( opc_SRCS
      oneParameterContinuation.cpp  )
SET( MY_EXECUTABLE
     "oneParameterContinuation.out" )
ADD_EXECUTABLE( ${MY_EXECUTABLE}
                ${opc_SRCS} )
SET_TARGET_PROPERTIES( ${MY_EXECUTABLE}
                       PROPERTIES LINKER_LANGUAGE Fortran )
TARGET_LINK_LIBRARIES( ${MY_EXECUTABLE}
                       ${internal_LIBS}
                       ${Trilinos_LIBRARIES}
                       ${Trilinos_TPL_LIBRARIES}
                       ${Boost_LIBRARIES} )

#SET ( opcConstraintInterface_SRCS
#      opcConstraintInterface.cpp  )
#SET( MY_EXECUTABLE
#     "opcConstraintInterface.out" )
#ADD_EXECUTABLE( ${MY_EXECUTABLE}
#                ${opcConstraintInterface_SRCS} )
#SET_TARGET_PROPERTIES( ${MY_EXECUTABLE}
#                       PROPERTIES LINKER_LANGUAGE Fortran )
#TARGET_LINK_LIBRARIES( ${MY_EXECUTABLE}
#                       ${internal_LIBS}
#                       ${Trilinos_LIBRARIES}
#                       ${Trilinos_TPL_LIBRARIES}
#                       ${Boost_LIBRARIES} )
