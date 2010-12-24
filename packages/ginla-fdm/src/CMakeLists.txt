SET ( ginlafdm_SRCS
      Ginla_FDM_State.cpp
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
      Ginla_State_Virtual.cpp
      Ginla_CreateSavable_Virtual.cpp
      Ginla_State_Updatable.cpp
      Ginla_StateTranslator_Virtual.cpp
      Ginla_ParameterHost_Virtual.cpp
      Ginla_Helpers.cpp
      Ginla_MagneticVectorPotential_Virtual.cpp
      Ginla_MagneticVectorPotential_X.cpp
      Ginla_MagneticVectorPotential_Y.cpp
      Ginla_MagneticVectorPotential_Z.cpp
      Ginla_MagneticVectorPotential_MagneticDot.cpp
      Ginla_MagneticVectorPotential_ZSquareSymmetric.cpp
      Ginla_IO_SaveNewtonData.cpp
      Ginla_IO_SaveEigenData.cpp
      Ginla_IO_StatsWriter.cpp
      Ginla_IO_StateWriter.cpp
      Ginla_IO_NoxObserver.cpp
      Ginla_StatusTest_Energy.cpp
      Ginla_StatusTest_Loop.cpp
      Ginla_StatusTest_Turnaround.cpp
      Ginla_StatusTest_MaxAcceptedSteps.cpp
      Ginla_StatusTest_ParameterLimits.cpp
      Ginla_StatusTest_StabilityChange.cpp
      Ginla_StatusTest_Factory.cpp
    )

FIND_PACKAGE( Trilinos REQUIRED )
INCLUDE_DIRECTORIES( ${Trilinos_INCLUDE_DIRS} )

INCLUDE_DIRECTORIES( ${ginzburg-landau_SOURCE_DIR}/packages/ginla-fdm/src
                     ${ginzburg-landau_SOURCE_DIR}/packages/ginla-fdm/src-komplex2
                     ${ginzburg-landau_SOURCE_DIR}/packages/ginla-fdm/src-recti
                     ${ginzburg-landau_SOURCE_DIR}/packages/ginla-fdm/src-vio
                   )

# workaround for Teuchos not appending Boost INCLUDES
FIND_PACKAGE( Boost REQUIRED )
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIRS} )

#ADD_LIBRARY ( ginlafdm STATIC
ADD_LIBRARY ( ginlafdm
              ${ginlafdm_SRCS} )

TARGET_LINK_LIBRARIES ( ginlafdm )
