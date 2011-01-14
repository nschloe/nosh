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
SET ( keoaztecoo_SRCS
      keo-aztecoo.cpp  )
SET ( keoaztecoo_LIBS
      ginlafvm
    )
SET( MY_EXECUTABLE
     "keo-aztecoo.exe" )
ADD_EXECUTABLE( ${MY_EXECUTABLE}
                ${keoaztecoo_SRCS} )
#SET_TARGET_PROPERTIES( ${MY_EXECUTABLE}
#                       PROPERTIES LINKER_LANGUAGE Fortran )
TARGET_LINK_LIBRARIES( ${MY_EXECUTABLE}
                       ${keoaztecoo_LIBS}
                       ${Trilinos_LIBRARIES}
                       ${Trilinos_TPL_LIBRARIES}
                       ${Boost_LIBRARIES} )
# ------------------------------------------------------------------------------
SET ( keobelos_SRCS
      keo-belos.cpp  )
SET ( keobelos_LIBS
      ginlafvm
    )
SET( MY_EXECUTABLE
     "keo-belos.exe" )
ADD_EXECUTABLE( ${MY_EXECUTABLE}
                ${keobelos_SRCS} )
#SET_TARGET_PROPERTIES( ${MY_EXECUTABLE}
#                       PROPERTIES LINKER_LANGUAGE Fortran )
TARGET_LINK_LIBRARIES( ${MY_EXECUTABLE}
                       ${keobelos_LIBS}
                       ${Trilinos_LIBRARIES}
                       ${Trilinos_TPL_LIBRARIES}
                       gfortran
                       ${Boost_LIBRARIES} )
# ------------------------------------------------------------------------------
SET ( jacobiansolver_SRCS
      jacobian-solver.cpp  )
SET ( jacobiansolver_LIBS
      ginlafvm
    )
SET( MY_EXECUTABLE
     "jacobian-solver.exe" )
ADD_EXECUTABLE( ${MY_EXECUTABLE}
                ${jacobiansolver_SRCS} )
#SET_TARGET_PROPERTIES( ${MY_EXECUTABLE}
#                       PROPERTIES LINKER_LANGUAGE Fortran )
TARGET_LINK_LIBRARIES( ${MY_EXECUTABLE}
                       ${jacobiansolver_LIBS}
                       ${Trilinos_LIBRARIES}
                       ${Trilinos_TPL_LIBRARIES}
                       gfortran
                       ${Boost_LIBRARIES} )
# ------------------------------------------------------------------------------
SET ( matvec_SRCS
      keo-matvec.cpp  )
SET ( matvec_LIBS
      ginlafvm
    )
SET( MY_EXECUTABLE
     "keo-matvec.exe" )
ADD_EXECUTABLE( ${MY_EXECUTABLE}
                ${matvec_SRCS} )
#SET_TARGET_PROPERTIES( ${MY_EXECUTABLE}
#                       PROPERTIES LINKER_LANGUAGE Fortran )
TARGET_LINK_LIBRARIES( ${MY_EXECUTABLE}
                       ${matvec_LIBS}
                       ${Trilinos_LIBRARIES}
                       ${Trilinos_TPL_LIBRARIES}
                       ${Boost_LIBRARIES} )
# ------------------------------------------------------------------------------
