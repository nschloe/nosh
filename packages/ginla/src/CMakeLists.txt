SET ( ginla_SRCS
      Ginla_State_Virtual.cpp
      Ginla_ParameterHost_Virtual.cpp
      Ginla_Helpers.cpp
      Ginla_MagneticVectorPotential_Centered.cpp
      Ginla_IO_SaveNewtonData.cpp
      Ginla_IO_SaveEigenData.cpp
      Ginla_IO_StatsWriter.cpp
      Ginla_IO_StateWriter.cpp
      Ginla_IO_NoxObserver.cpp
      Ginla_Komplex_DoubleMatrix.cpp
      Ginla_Komplex_LinearProblem.cpp
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

#INCLUDE_DIRECTORIES( ${ginzburg-landau_SOURCE_DIR}/packages/recti/src
#                     ${ginzburg-landau_SOURCE_DIR}/packages/vio/src )

#ADD_LIBRARY ( ginla STATIC
ADD_LIBRARY ( ginla
              ${ginla_SRCS} )
