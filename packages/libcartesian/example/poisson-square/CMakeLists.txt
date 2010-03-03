CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

FIND_PACKAGE( NOX  REQUIRED )
INCLUDE_DIRECTORIES ( ${NOX_INCLUDE_DIRS} )
LINK_DIRECTORIES( ${Trilinos_LIBRARY_DIRS}
                  ${Trilinos_TPL_LIBRARY_DIRS} )

FIND_PACKAGE( Boost COMPONENTS filesystem REQUIRED )
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIRS} )
SET(Boost_USE_STATIC_LIBS   ON)
SET(Boost_USE_MULTITHREADED ON)
LINK_DIRECTORIES( ${Boost_LIBRARY_DIRS} )


INCLUDE_DIRECTORIES ( ${ginzburg-landau_SOURCE_DIR}/packages/libio/src
                      ${ginzburg-landau_SOURCE_DIR}/packages/libcartesian/src
                    )

SET ( poisson_square_SRCS
      main.cpp )

SET ( internal_LIBS
      cartesian
      io
    )

SET ( MY_EXECUTABLE "poisson-square.out" )

ADD_EXECUTABLE ( ${MY_EXECUTABLE}
                 ${poisson_square_SRCS} )

TARGET_LINK_LIBRARIES ( ${MY_EXECUTABLE}
                        ${internal_LIBS}
                       ${Trilinos_LIBRARIES}
                       ${Trilinos_TPL_LIBRARIES}
                       ${Boost_LIBRARIES}
                       ${GFORTRAN} )
                      )