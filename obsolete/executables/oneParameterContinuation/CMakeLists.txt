CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

FIND_PACKAGE( NOX REQUIRED )
LINK_DIRECTORIES( ${Trilinos_LIBRARY_DIRS}
                  ${Trilinos_TPL_LIBRARY_DIRS} )
INCLUDE_DIRECTORIES ( ${Trilinos_INCLUDE_DIRS}
                      ${Trilinos_TPL_INCLUDE_DIRS} )
                      
FIND_PACKAGE( Boost COMPONENTS system filesystem REQUIRED )
INCLUDE_DIRECTORIES ( ${Boost_INCLUDE_DIRS} )
LINK_DIRECTORIES( ${Boost_LIBRARY_DIRS} )


SET ( internal_LIBS
      ginlafdm
      ginla
      recti
      vio
    )

INCLUDE_DIRECTORIES ( ${ginzburg-landau_SOURCE_DIR}/packages/vio/src
                      ${ginzburg-landau_SOURCE_DIR}/packages/recti/src
                      ${ginzburg-landau_SOURCE_DIR}/packages/ginla/src
                      ${ginzburg-landau_SOURCE_DIR}/packages/ginla/src-fdm
                    )

SET ( opc_SRCS
      oneParameterContinuation.cpp  )
SET( MY_EXECUTABLE
     "oneParameterContinuation.out" )
ADD_EXECUTABLE( ${MY_EXECUTABLE}
                ${opc_SRCS} )
TARGET_LINK_LIBRARIES( ${MY_EXECUTABLE}
                       ${internal_LIBS}
                       ${Trilinos_LIBRARIES}
                       ${Trilinos_TPL_LIBRARIES}
                       ${Boost_LIBRARIES} )

SET ( opc_SRCS
      opcWithModelEvaluator.cpp  )
SET( MY_EXECUTABLE
     "opcWithModelEvaluator.out" )
ADD_EXECUTABLE( ${MY_EXECUTABLE}
                ${opc_SRCS} )
TARGET_LINK_LIBRARIES( ${MY_EXECUTABLE}
                       ${internal_LIBS}
                       ${Trilinos_LIBRARIES}
                       ${Trilinos_TPL_LIBRARIES}
                       ${Boost_LIBRARIES} )

# SET ( opcConstraintInterface_SRCS
#       opcConstraintInterface.cpp  )
# SET( MY_EXECUTABLE
#      "opcConstraintInterface.out" )
# ADD_EXECUTABLE( ${MY_EXECUTABLE}
#                 ${opcConstraintInterface_SRCS} )
# TARGET_LINK_LIBRARIES( ${MY_EXECUTABLE}
#                        ${internal_LIBS}
#                        ${Trilinos_LIBRARIES}
#                        ${Trilinos_TPL_LIBRARIES}
#                        ${Boost_LIBRARIES} )
# 
# SET ( opcModelEvaluator_SRCS
#       opcModelEvaluator.cpp  )
# SET( MY_EXECUTABLE
#      "opcModelEvaluator.out" )
# ADD_EXECUTABLE( ${MY_EXECUTABLE}
#                 ${opcModelEvaluator_SRCS} )
# TARGET_LINK_LIBRARIES( ${MY_EXECUTABLE}
#                        ${internal_LIBS}
#                        ${Trilinos_LIBRARIES}
#                        ${Trilinos_TPL_LIBRARIES}
#                        ${Boost_LIBRARIES} )
