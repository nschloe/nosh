SET ( ginla_SRCS
      Ginla_EpetraFVM_State.cpp
      Ginla_EpetraFVM_ModelEvaluator.cpp
      Ginla_EpetraFVM_JacobianOperator.cpp
      Ginla_EpetraFVM_KineticEnergyOperator.cpp
      Ginla_EpetraFVM_KeoPreconditioner.cpp
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
