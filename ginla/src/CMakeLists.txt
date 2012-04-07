SET ( ginla_SRCS
      dummy.f
      Ginla_Helpers.cpp
      Ginla_State.cpp
      Ginla_ModelEvaluator.cpp
      Ginla_JacobianOperator.cpp
      Ginla_KeoContainer.cpp
      Ginla_KeoRegularized.cpp
      Ginla_MagneticVectorPotential_Virtual.cpp
      Ginla_MagneticVectorPotential_ExplicitValues.cpp
      Ginla_MagneticVectorPotential_ConstantInSpace.cpp
      Ginla_StkMesh.cpp
      Ginla_StkMeshReader.cpp
      Ginla_StatsWriter.cpp
      Ginla_NoxObserver.cpp
      Ginla_SaveEigenData.cpp
    )

FIND_PACKAGE( Trilinos REQUIRED )
INCLUDE_DIRECTORIES( ${Trilinos_INCLUDE_DIRS}
                     ${Trilinos_TPL_INCLUDE_DIRS} )

ADD_LIBRARY ( ginla
              ${ginla_SRCS} )
