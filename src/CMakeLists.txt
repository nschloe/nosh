SET ( cuantico_SRCS
      dummy.f
      Cuantico_Helpers.cpp
      Cuantico_State.cpp
      Cuantico_ModelEvaluator.cpp
      Cuantico_JacobianOperator.cpp
      Cuantico_KeoContainer.cpp
      Cuantico_KeoRegularized.cpp
      Cuantico_ScalarPotential_Virtual.cpp
      Cuantico_ScalarPotential_Constant.cpp
      Cuantico_MagneticVectorPotential_Virtual.cpp
      Cuantico_MagneticVectorPotential_ExplicitValues.cpp
      Cuantico_MagneticVectorPotential_ConstantField.cpp
      Cuantico_StkMesh.cpp
      Cuantico_StkMeshReader.cpp
      Cuantico_CsvWriter.cpp
      Cuantico_NoxObserver.cpp
      Cuantico_SaveEigenData.cpp
    )

FIND_PACKAGE( Trilinos REQUIRED )
INCLUDE_DIRECTORIES( ${Trilinos_INCLUDE_DIRS}
                     ${Trilinos_TPL_INCLUDE_DIRS} )

ADD_LIBRARY(cuantico
            ${cuantico_SRCS} )
