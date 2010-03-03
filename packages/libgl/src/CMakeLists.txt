ENABLE_LANGUAGE( CXX )

SET ( gl_SRCS
      GL_MagneticVectorPotential_Centered.cpp
      ginzburgLandau.cpp
      GL_Operator_Virtual.cpp
      GL_Operator_BCInner.cpp
      GL_Operator_BCOuter.cpp
      GL_Operator_BCCentral.cpp
      GL_IO_SaveNewtonData.cpp
      GL_IO_SaveEigenData.cpp
      GL_Komplex.cpp
      GL_LocaSystem_Default.cpp
      GL_LocaSystem_Bordered.cpp
      GL_LocaSystem_Predictor.cpp
      AbstractStateWriter.cpp
    )

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
