# actually only needed: NOX and Tpetra
FIND_PACKAGE( Trilinos REQUIRED )
INCLUDE_DIRECTORIES ( ${Trilinos_INCLUDE_DIRS} )
LINK_DIRECTORIES( ${Trilinos_LIBRARY_DIRS} )

FIND_PACKAGE( Boost COMPONENTS filesystem REQUIRED )
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIRS} )
# SET(Boost_USE_STATIC_LIBS   ON)
SET(Boost_USE_MULTITHREADED ON)

INCLUDE_DIRECTORIES ( ${ginzburg-landau_SOURCE_DIR}/packages/vio/src
                      ${ginzburg-landau_SOURCE_DIR}/packages/recti/src
                    )

SET ( poisson_SRCS
      main.cpp )

SET( MY_EXECUTABLE
     "poisson.out" )

ADD_EXECUTABLE ( ${MY_EXECUTABLE}
                 ${poisson_SRCS} )

SET_TARGET_PROPERTIES( ${MY_EXECUTABLE}
                       PROPERTIES LINKER_LANGUAGE Fortran )

SET ( internal_LIBS
      recti
      vio
    )

# Set executable linking information.
TARGET_LINK_LIBRARIES ( ${MY_EXECUTABLE}
                        ${internal_LIBS}
                        ${Trilinos_LIBRARIES}
                        ${Trilinos_TPL_LIBRARIES}
                        ${Boost_LIBRARIES}
                        ${LINK_LIBRARIES}
                      )
