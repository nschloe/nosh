SET ( ginla_SRCS
      Ginla_Helpers.cpp
      Ginla_MagneticVectorPotential_Centered.cpp
      Ginla_Operator_Virtual.cpp
      Ginla_Operator_BCInner.cpp
      Ginla_Operator_BCOuter.cpp
      Ginla_Operator_BCCentral.cpp
      Ginla_IO_SaveNewtonData.cpp
      Ginla_IO_SaveEigenData.cpp
      Ginla_IO_StatsWriter.cpp
      Ginla_IO_StateWriter.cpp
      Ginla_Komplex.cpp
      Ginla_LocaSystem_Default.cpp
      Ginla_LocaSystem_Bordered.cpp
      Ginla_Constraint_MinDist.cpp
#       Ginla_Constraint_MinDistReal.cpp
      Ginla_Perturbation_Quadrants.cpp
      Ginla_Perturbation_Virtual.cpp
    )

FIND_PACKAGE( Trilinos REQUIRED )
INCLUDE_DIRECTORIES( ${Trilinos_INCLUDE_DIRS} )

INCLUDE_DIRECTORIES( ${ginzburg-landau_SOURCE_DIR}/packages/recti/src
                     ${ginzburg-landau_SOURCE_DIR}/packages/vio/src )

ADD_LIBRARY ( ginla STATIC
              ${ginla_SRCS} )
