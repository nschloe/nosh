SET ( ginla_SRCS
      Ginla_FVM_State.cpp
      Ginla_FVM_ModelEvaluator.cpp
    )

FIND_PACKAGE( Trilinos REQUIRED )
INCLUDE_DIRECTORIES( ${Trilinos_INCLUDE_DIRS} )

FIND_PACKAGE( VTK REQUIRED )
INCLUDE_DIRECTORIES( ${VTK_INCLUDE_DIRS} )

INCLUDE_DIRECTORIES( ${ginzburg-landau_SOURCE_DIR}/packages/ginla-fvm2/src-vio2
                     ${ginzburg-landau_SOURCE_DIR}/packages/ginla-fdm/src-komplex2
                   )

ADD_LIBRARY ( ginlafvm-tpetra
              ${ginla_SRCS} )
