SET ( komplex2_SRCS
      Komplex2_DoubleMatrix.cpp
      Komplex2_LinearProblem.cpp
   )

FIND_PACKAGE( Trilinos REQUIRED )
INCLUDE_DIRECTORIES( ${Trilinos_INCLUDE_DIRS} )

ADD_LIBRARY ( komplex2
              ${komplex2_SRCS} )
