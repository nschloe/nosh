SET ( ginla_SRCS
      dummy.f
      Ginla_Helpers.cpp
      Ginla_EpetraFVM_State.cpp
      Ginla_EpetraFVM_ModelEvaluator.cpp
      Ginla_EpetraFVM_JacobianOperator.cpp
      Ginla_EpetraFVM_KeoFactory.cpp
      Ginla_EpetraFVM_KeoPreconditioner.cpp
      Ginla_MagneticVectorPotential.cpp
      Ginla_EpetraFVM_StkMesh.cpp
      Ginla_EpetraFVM_StkMeshReader.cpp
      Ginla_EpetraFVM_StkMeshWriter.cpp
      Ginla_IO_StateWriter.cpp
      Ginla_IO_StatsWriter.cpp
      Ginla_IO_NoxObserver.cpp
      Ginla_IO_SaveEigenData.cpp
    )

FIND_PACKAGE( Trilinos REQUIRED )
INCLUDE_DIRECTORIES( ${Trilinos_INCLUDE_DIRS} )

ADD_LIBRARY ( ginla
              ${ginla_SRCS} )
