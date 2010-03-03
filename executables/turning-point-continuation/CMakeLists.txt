CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

FIND_PACKAGE( NOX REQUIRED PATHS /opt/trilinos/dev/gcc/4.3.4/ )
FIND_PACKAGE( VTK REQUIRED )
FIND_PACKAGE( Boost COMPONENTS system filesystem REQUIRED )
SET(Boost_USE_STATIC_LIBS   ON)
SET(Boost_USE_MULTITHREADED ON)

FIND_LIBRARY( GFORTRAN ${GFORTRAN_LIBRARY_NAMES}
              HINTS ${GFORTRAN_LIBRARY_DIRS}    )

SET ( ginzburg_landau_SRCS
      glTurningPointContinuation.cpp )

SET ( internal_LIBS
      gl
      cartesian
      io
    )

INCLUDE_DIRECTORIES ( ${ginzburg-landau_SOURCE_DIR}/packages/libio/src
                      ${ginzburg-landau_SOURCE_DIR}/packages/libcartesian/src
                      ${ginzburg-landau_SOURCE_DIR}/packages/libgl/src
                      ${VTK_INCLUDE_DIRS}
                      ${Trilinos_INCLUDE_DIRS}
                      ${Trilinos_TPL_INCLUDE_DIRS} )

LINK_DIRECTORIES( ${Boost_LIBRARY_DIRS}
                  ${VTK_LIBRARY_DIRS}
                  ${Trilinos_LIBRARY_DIRS}
                  ${Trilinos_TPL_LIBRARY_DIRS} )

SET( MY_EXECUTABLE "glTurningPointContinuation.out" )

ADD_EXECUTABLE( ${MY_EXECUTABLE}
                ${ginzburg_landau_SRCS})

TARGET_LINK_LIBRARIES( ${MY_EXECUTABLE}
                       ${internal_LIBS}
                       ${Trilinos_LIBRARIES}
                       ${Trilinos_TPL_LIBRARIES}
                       ${Boost_LIBRARIES}
                       ${GFORTRAN} )