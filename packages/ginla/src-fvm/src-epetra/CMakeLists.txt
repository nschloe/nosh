SET ( ginla_SRCS
      Ginla_FVM_State.cpp
      Ginla_FVM_ModelEvaluator.cpp
      Ginla_FVM_JacobianOperator.cpp
      Ginla_FVM_KineticEnergy.cpp
    )

FIND_PACKAGE( Trilinos REQUIRED )
INCLUDE_DIRECTORIES( ${Trilinos_INCLUDE_DIRS} )

FIND_PACKAGE( VTK REQUIRED )
INCLUDE_DIRECTORIES( ${VTK_INCLUDE_DIRS} )

INCLUDE_DIRECTORIES( ${ginzburg-landau_SOURCE_DIR}/packages/ginla/src
                     ${ginzburg-landau_SOURCE_DIR}/packages/komplex2/src
                     ${ginzburg-landau_SOURCE_DIR}/packages/vio/src
                   )

ADD_LIBRARY ( ginlafvm
              ${ginla_SRCS} )
