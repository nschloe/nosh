CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

INCLUDE_DIRECTORIES ( ${ginzburg-landau_SOURCE_DIR}/src
                      ${Trilinos_INCLUDE_DIRS}
                      ${Trilinos_TPL_INCLUDE_DIRS} )

LINK_DIRECTORIES( ${Trilinos_LIBRARY_DIRS}
                  ${Trilinos_TPL_LIBRARY_DIRS} )
# ------------------------------------------------------------------------------
SET ( keoaztecoo_SRCS
      keo-aztecoo.cpp  )
SET ( keoaztecoo_LIBS
      ginla
    )
SET( MY_EXECUTABLE
     "keo-aztecoo.exe" )
ADD_EXECUTABLE( ${MY_EXECUTABLE}
                ${keoaztecoo_SRCS} )
TARGET_LINK_LIBRARIES( ${MY_EXECUTABLE}
                       ${keoaztecoo_LIBS}
                       ${Trilinos_LIBRARIES}
                       ${Trilinos_TPL_LIBRARIES} )
# ------------------------------------------------------------------------------
SET ( keobelos_SRCS
      keo-belos.cpp  )
SET ( keobelos_LIBS
      ginla
    )
SET( MY_EXECUTABLE
     "keo-belos.exe" )
ADD_EXECUTABLE( ${MY_EXECUTABLE}
                ${keobelos_SRCS} )
TARGET_LINK_LIBRARIES( ${MY_EXECUTABLE}
                       ${keobelos_LIBS}
                       ${Trilinos_LIBRARIES}
                       ${Trilinos_TPL_LIBRARIES} )
# ------------------------------------------------------------------------------
SET ( jacobiansolver_SRCS
      jacobian-solver.cpp  )
SET ( jacobiansolver_LIBS
      ginla
    )
SET( MY_EXECUTABLE
     "jacobian-solver.exe" )
ADD_EXECUTABLE( ${MY_EXECUTABLE}
                ${jacobiansolver_SRCS} )
TARGET_LINK_LIBRARIES( ${MY_EXECUTABLE}
                       ${jacobiansolver_LIBS}
                       ${Trilinos_LIBRARIES}
                       ${Trilinos_TPL_LIBRARIES} )
# ------------------------------------------------------------------------------
SET ( matvec_SRCS
      keo-matvec.cpp  )
SET ( matvec_LIBS
      ginla
    )
SET( MY_EXECUTABLE
     "keo-matvec.exe" )
ADD_EXECUTABLE( ${MY_EXECUTABLE}
                ${matvec_SRCS} )
TARGET_LINK_LIBRARIES( ${MY_EXECUTABLE}
                       ${matvec_LIBS}
                       ${Trilinos_LIBRARIES}
                       ${Trilinos_TPL_LIBRARIES} )
# ------------------------------------------------------------------------------
