INCLUDE_DIRECTORIES(${Nosh_SOURCE_DIR}/src
                    ${Trilinos_INCLUDE_DIRS}
                    ${Trilinos_TPL_INCLUDE_DIRS})

LINK_DIRECTORIES( ${Boost_LIBRARY_DIRS}
                  ${Trilinos_LIBRARY_DIRS}
                  ${Trilinos_TPL_LIBRARY_DIRS} )
# ------------------------------------------------------------------------------
SET(MY_EXECUTABLE "noshExamples")
ADD_EXECUTABLE( ${MY_EXECUTABLE} "exampleCalls.cpp")
TARGET_LINK_LIBRARIES(${MY_EXECUTABLE}
                      "nosh"
                      ${Trilinos_LIBRARIES}
                      ${Trilinos_TPL_LIBRARIES})

# Don't install this example code.
# ------------------------------------------------------------------------------
