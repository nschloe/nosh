INCLUDE( ${ginzburg-landau_SOURCE_DIR}/cmake-system.txt REQUIRED )

ENABLE_LANGUAGE( CXX )

SET ( gl_SRCS
      GridVirtual.cpp
      GridUniformVirtual.cpp
      GridSquare.cpp
      GridUniformSquare.cpp
      GridReader.cpp
      MagneticVectorPotential.cpp
      ginzburgLandau.cpp
      glBoundaryConditionsVirtual.cpp
      glBoundaryConditionsInner.cpp
      glBoundaryConditionsOuter.cpp
      glBoundaryConditionsCentral.cpp
      glPrePostOperator.cpp
      GlKomplex.cpp
      AbstractStateWriter.cpp
      glSystem.cpp
      glPredictorSystem.cpp
      GlSystemWithConstraint.cpp
      domainRectangle.cpp
      eigenSaver.cpp          )

INCLUDE( SetTrilinosVariables REQUIRED )
INCLUDE_DIRECTORIES( ${Trilinos_INCLUDE_DIR} )

IF(MPI_COMPILED)
    FIND_PACKAGE( MPI REQUIRED )
    INCLUDE_DIRECTORIES( ${MPI_INCLUDE_PATH} )
ENDIF(MPI_COMPILED)

INCLUDE_DIRECTORIES( ${ginzburg-landau_SOURCE_DIR}/packages/libio/src )

ADD_LIBRARY ( gl STATIC
              ${gl_SRCS} )

# list library dependencies here
#TARGET_LINK_LIBRARIES ( gl )
