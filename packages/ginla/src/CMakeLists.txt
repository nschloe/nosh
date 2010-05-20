SET ( ginla_SRCS
      Ginla_Helpers.cpp
      Ginla_MagneticVectorPotential_Centered.cpp
      Ginla_State.cpp
      Ginla_Operator_Virtual.cpp
#       Ginla_Operator_BCInner.cpp
#       Ginla_Operator_BCOuter.cpp
      Ginla_Operator_BCCentral.cpp
      Ginla_IO_SaveNewtonData.cpp
      Ginla_IO_SaveEigenData.cpp
      Ginla_IO_StatsWriter.cpp
      Ginla_IO_StateWriter.cpp
      Ginla_IO_NoxObserver.cpp
      Ginla_Komplex_DoubleMatrix.cpp
      Ginla_Komplex_LinearProblem.cpp
      Ginla_LocaSystem_Virtual.cpp
      Ginla_LocaSystem_Default.cpp
      Ginla_LocaSystem_Bordered.cpp
      Ginla_ModelEvaluator_Default.cpp
      Ginla_Constraint_MinDist.cpp
#       Ginla_Constraint_MinDistReal.cpp
      Ginla_Perturbation_Quadrants.cpp
      Ginla_Perturbation_Virtual.cpp
      Ginla_StatusTest_Energy.cpp
      Ginla_StatusTest_Loop.cpp
      Ginla_StatusTest_Turnaround.cpp
#      Ginla_BranchExplorer_Driver.cpp
#      Ginla_BranchExplorer_ProcessEigenData.cpp
    )

FIND_PACKAGE( Trilinos REQUIRED )
INCLUDE_DIRECTORIES( ${Trilinos_INCLUDE_DIRS} )

INCLUDE_DIRECTORIES( ${ginzburg-landau_SOURCE_DIR}/packages/recti/src
                     ${ginzburg-landau_SOURCE_DIR}/packages/vio/src )

#ADD_LIBRARY ( ginla STATIC
ADD_LIBRARY ( ginla
              ${ginla_SRCS} )
