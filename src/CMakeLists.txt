SET ( nosh_SRCS
      Nosh_Helpers.cpp
      Nosh_State.cpp
      Nosh_ModelEvaluator.cpp
      Nosh_JacobianOperator.cpp
      Nosh_KeoBuilder.cpp
      Nosh_KeoRegularized.cpp
      Nosh_ScalarPotential_Virtual.cpp
      Nosh_ScalarPotential_Constant.cpp
      Nosh_MagneticVectorPotential_Virtual.cpp
      Nosh_MagneticVectorPotential_ExplicitValues.cpp
      Nosh_MagneticVectorPotential_ConstantField.cpp
      Nosh_StkMesh.cpp
      Nosh_StkMeshReader.cpp
      Nosh_CsvWriter.cpp
      Nosh_NoxObserver.cpp
      Nosh_SaveEigenData.cpp
    )

FIND_PACKAGE( Trilinos REQUIRED )
INCLUDE_DIRECTORIES( ${Trilinos_INCLUDE_DIRS}
                     ${Trilinos_TPL_INCLUDE_DIRS} )

ADD_LIBRARY(nosh
            ${nosh_SRCS} )
