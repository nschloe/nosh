INCLUDE( ${ginzburg-landau_SOURCE_DIR}/cmake-system.txt REQUIRED )
ENABLE_LANGUAGE( CXX )

SET ( gl_SRCS
      MagneticVectorPotential.cpp
      ginzburgLandau.cpp
      GlOperatorVirtual.cpp
      GlOperatorBCInner.cpp
      GlOperatorBCOuter.cpp
      GlOperatorBCCentral.cpp
      glPrePostOperator.cpp
      GlKomplex.cpp
      GlHelpers.cpp
      AbstractStateWriter.cpp
      glSystem.cpp
      glPredictorSystem.cpp
      GlSystemWithConstraint.cpp
      eigenSaver.cpp          )

INCLUDE( SetTrilinosVariables REQUIRED )
INCLUDE_DIRECTORIES( ${Trilinos_INCLUDE_DIR} )

IF(MPI_COMPILED)
    FIND_PACKAGE( MPI REQUIRED )
    INCLUDE_DIRECTORIES( ${MPI_INCLUDE_PATH} )
ENDIF(MPI_COMPILED)

INCLUDE_DIRECTORIES( ${ginzburg-landau_SOURCE_DIR}/packages/libcartesian/src )
INCLUDE_DIRECTORIES( ${ginzburg-landau_SOURCE_DIR}/packages/libio/src )

ADD_LIBRARY ( gl STATIC
              ${gl_SRCS} )
