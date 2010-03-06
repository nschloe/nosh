SET ( gl_SRCS
      GL_Helpers.cpp
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
      GL_Constraint_MinDist.cpp
      GL_Perturbation_Quadrants.cpp
      GL_Perturbation_Virtual.cpp
      AbstractStateWriter.cpp
    )

FIND_PACKAGE( Trilinos REQUIRED )
INCLUDE_DIRECTORIES( ${Trilinos_INCLUDE_DIRS} )

INCLUDE_DIRECTORIES( ${ginzburg-landau_SOURCE_DIR}/packages/libcartesian/src
                     ${ginzburg-landau_SOURCE_DIR}/packages/libio/src )

ADD_LIBRARY ( gl STATIC
              ${gl_SRCS} )
