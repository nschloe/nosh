SET ( nosh_SRCS
      Nosh_Helpers.cpp
      Nosh_ModelEvaluator_Virtual.cpp
      Nosh_ModelEvaluator_Nls.cpp
      Nosh_ModelEvaluator_Bordered.cpp
      Nosh_JacobianOperator.cpp
      Nosh_MatrixBuilder_Virtual.cpp
      Nosh_MatrixBuilder_Keo.cpp
      Nosh_MatrixBuilder_Laplace.cpp
      Nosh_KeoRegularized.cpp
      Nosh_ScalarField_Virtual.cpp
      Nosh_ScalarField_Constant.cpp
      Nosh_ScalarField_ConstantWithParam.cpp
      Nosh_VectorField_Virtual.cpp
      Nosh_VectorField_ExplicitValues.cpp
      Nosh_VectorField_ConstantCurl.cpp
      Nosh_StkMesh.cpp
      Nosh_StkMeshReader.cpp
      Nosh_CsvWriter.cpp
      Nosh_Observer.cpp
      Nosh_SaveEigenData.cpp
    )

FIND_PACKAGE( Trilinos REQUIRED )
INCLUDE_DIRECTORIES( ${Trilinos_INCLUDE_DIRS}
                     ${Trilinos_TPL_INCLUDE_DIRS} )

ADD_LIBRARY(nosh
            ${nosh_SRCS} )
