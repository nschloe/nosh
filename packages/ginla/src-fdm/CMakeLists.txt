SET ( ginlafdm_SRCS
      Ginla_FDM_State.cpp
      Ginla_FDM_Constraint_MinDist.cpp
#       Ginla_FDM_Constraint_MinDistReal.cpp
      Ginla_FDM_LocaSystem_Virtual.cpp
      Ginla_FDM_LocaSystem_Default.cpp
      Ginla_FDM_LocaSystem_Bordered.cpp
      Ginla_FDM_ModelEvaluator_Default.cpp
      Ginla_FDM_ModelEvaluator_Bordered.cpp
      Ginla_FDM_Operator_Virtual.cpp
      Ginla_FDM_Operator_BCCentral.cpp
#       Ginla_FDM_Operator_BCInner.cpp
#       Ginla_FDM_Operator_BCOuter.cpp
#       Ginla_FDM_Perturbation_Virtual.cpp
#       Ginla_FDM_Perturbation_Quadrants.cpp
#       Ginla_FDM_Preconditioner_Virtual.cpp
#       Ginla_FDM_Preconditioner_PerturbedNeumannLaplace.cpp
    )

FIND_PACKAGE( Trilinos REQUIRED )
INCLUDE_DIRECTORIES( ${Trilinos_INCLUDE_DIRS} )

INCLUDE_DIRECTORIES( ${ginzburg-landau_SOURCE_DIR}/packages/ginla/src
                     ${ginzburg-landau_SOURCE_DIR}/packages/recti/src
                     ${ginzburg-landau_SOURCE_DIR}/packages/vio/src
                   )

#ADD_LIBRARY ( ginlafdm STATIC
ADD_LIBRARY ( ginlafdm
              ${ginlafdm_SRCS} )

TARGET_LINK_LIBRARIES ( ginlafdm
                        ginla )