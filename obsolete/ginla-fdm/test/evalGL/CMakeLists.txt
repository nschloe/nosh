FIND_PACKAGE( Trilinos REQUIRED )
INCLUDE_DIRECTORIES ( ${Trilinos_INCLUDE_DIRS}
                      ${Trilinos_TPL_INCLUDE_DIRS} )
LINK_DIRECTORIES( ${Trilinos_LIBRARY_DIRS}
                  ${Trilinos_TPL_LIBRARY_DIRS} )

FIND_PACKAGE( Boost
              COMPONENTS system filesystem unit_test_framework
              REQUIRED )
INCLUDE_DIRECTORIES ( ${Boost_INCLUDE_DIRS} )

INCLUDE_DIRECTORIES ( ${ginzburg-landau_SOURCE_DIR}/packages/vio/src
                      ${ginzburg-landau_SOURCE_DIR}/packages/recti/src
                      ${ginzburg-landau_SOURCE_DIR}/packages/ginla/src
                      ${ginzburg-landau_SOURCE_DIR}/packages/ginla/src-fdm
                    )

SET( MY_EXECUTABLE
     "evalGl.out" )

# see http://www.docstoc.com/docs/14857835/CTest-and-Dashboard-_Slides-by-J
#INCLUDE( CTest )
ADD_TEST( evalGl0
          ${MY_EXECUTABLE}
          "--input-state-file=${CMAKE_CURRENT_SOURCE_DIR}/data/input0.vti"
          "--reference-residual-file=${CMAKE_CURRENT_SOURCE_DIR}/data/reference0.vti"
          "--reference-jacobian-file=${CMAKE_CURRENT_SOURCE_DIR}/data/jacobian0.h5"
        )

ADD_TEST( evalGl1
          ${MY_EXECUTABLE}
          "--input-state-file=${CMAKE_CURRENT_SOURCE_DIR}/data/input1.vti"
          "--reference-residual-file=${CMAKE_CURRENT_SOURCE_DIR}/data/reference1.vti"
          "--reference-jacobian-file=${CMAKE_CURRENT_SOURCE_DIR}/data/jacobian1.h5"
        )

ADD_TEST( evalGl2
          ${MY_EXECUTABLE}
          "--input-state-file=${CMAKE_CURRENT_SOURCE_DIR}/data/input2.vti"
          "--reference-residual-file=${CMAKE_CURRENT_SOURCE_DIR}/data/reference2.vti"
          "--reference-jacobian-file=${CMAKE_CURRENT_SOURCE_DIR}/data/jacobian2.h5"
        )

# ADD_TEST( zeroStepLocaMultiproc
#           "/usr/bin/env" "mpirun" "-n" "2"
#           "./${MY_EXECUTABLE}"
#           "--xml-input-file=${CMAKE_CURRENT_SOURCE_DIR}/data/conf.xml"
#           "--expected-solution-file=${CMAKE_CURRENT_SOURCE_DIR}/data/expected-solution.vti"
#         )

ADD_EXECUTABLE ( ${MY_EXECUTABLE}
                 main.cpp )

SET ( internal_LIBS
      ginlafdm
      recti
      vio
    )

# Set executable linking information.
TARGET_LINK_LIBRARIES ( ${MY_EXECUTABLE}
                        ${internal_LIBS}
                        ${Trilinos_LIBRARIES}
                        ${Trilinos_TPL_LIBRARIES}
                        ${Boost_LIBRARIES}
                      )
