CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

FIND_PACKAGE( Trilinos REQUIRED )
INCLUDE_DIRECTORIES ( ${Trilinos_INCLUDE_DIRS}
                      ${Trilinos_TPL_INCLUDE_DIRS} )
LINK_DIRECTORIES( ${Trilinos_LIBRARY_DIRS}
                  ${Trilinos_TPL_LIBRARY_DIRS} )

FIND_PACKAGE( VTK REQUIRED )

FIND_PACKAGE( Boost COMPONENTS system filesystem REQUIRED )
INCLUDE_DIRECTORIES ( ${Boost_INCLUDE_DIRS} )

SET ( internal_LIBS
      ginla
      recti
      vio
    )

INCLUDE_DIRECTORIES ( ${ginzburg-landau_SOURCE_DIR}/packages/vio/src
                      ${ginzburg-landau_SOURCE_DIR}/packages/recti/src
                      ${ginzburg-landau_SOURCE_DIR}/packages/ginla/src
                      ${VTK_INCLUDE_DIRS} )

LINK_DIRECTORIES( ${Boost_LIBRARY_DIRS}
                  ${VTK_LIBRARY_DIRS} )

SET( MY_EXECUTABLE
     "gl-nox.out" )
SET ( glNOX_SRCS
      glNoxHelpers.cpp
      main.cpp )
ADD_EXECUTABLE( ${MY_EXECUTABLE}
                ${glNOX_SRCS})

SET_TARGET_PROPERTIES( ${MY_EXECUTABLE}
                       PROPERTIES LINKER_LANGUAGE Fortran )

TARGET_LINK_LIBRARIES( ${MY_EXECUTABLE}
                       ${internal_LIBS}
                       ${Trilinos_LIBRARIES}
                       ${Trilinos_TPL_LIBRARIES}
                       ${Boost_LIBRARIES} )

SET( MY_EXECUTABLE
     "noxWithModelEvaluator.out" )
SET ( glNOXModelEvaluator_SRCS
      noxWithModelEvaluator.cpp )
ADD_EXECUTABLE( ${MY_EXECUTABLE}
                ${glNOXModelEvaluator_SRCS})
SET_TARGET_PROPERTIES( ${MY_EXECUTABLE}
                       PROPERTIES LINKER_LANGUAGE Fortran )
TARGET_LINK_LIBRARIES( ${MY_EXECUTABLE}
                       ${internal_LIBS}
                       ${Trilinos_LIBRARIES}
                       ${Trilinos_TPL_LIBRARIES}
                       ${Boost_LIBRARIES} )
