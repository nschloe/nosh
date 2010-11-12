SET ( ginla_SRCS
      Ginla_Helpers.cpp
      Ginla_EpetraFVM_State.cpp
      Ginla_EpetraFVM_ModelEvaluator.cpp
      Ginla_EpetraFVM_JacobianOperator.cpp
      Ginla_EpetraFVM_KeoFactory.cpp
      Ginla_EpetraFVM_KeoPreconditioner.cpp
      Ginla_EpetraFVM_StkMesh.cpp
      Ginla_EpetraFVM_StkMeshReader.cpp
      Ginla_EpetraFVM_StkMeshWriter.cpp
      Ginla_MagneticVectorPotential_Virtual.cpp
      Ginla_MagneticVectorPotential_X.cpp
      Ginla_MagneticVectorPotential_Y.cpp
      Ginla_MagneticVectorPotential_Z.cpp
      Ginla_MagneticVectorPotential_ZSquareSymmetric.cpp
      Ginla_MagneticVectorPotential_MagneticDot.cpp
      Ginla_IO_StateWriter.cpp
      Ginla_IO_StatsWriter.cpp
      Ginla_IO_NoxObserver.cpp
      Ginla_IO_SaveEigenData.cpp
    )

FIND_PACKAGE( Trilinos REQUIRED )
INCLUDE_DIRECTORIES( ${Trilinos_INCLUDE_DIRS} )

FIND_PACKAGE( VTK REQUIRED )
INCLUDE_DIRECTORIES( ${VTK_INCLUDE_DIRS} )

INCLUDE_DIRECTORIES( ${ginzburg-landau_SOURCE_DIR}/packages/ginla/src
                     ${ginzburg-landau_SOURCE_DIR}/packages/komplex2/src
                     ${ginzburg-landau_SOURCE_DIR}/packages/vio/src
                   )

ADD_LIBRARY ( ginlafvm-epetra
              ${ginla_SRCS} )
