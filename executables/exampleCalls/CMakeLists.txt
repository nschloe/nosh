INCLUDE_DIRECTORIES(${Nosh_SOURCE_DIR})
INCLUDE_DIRECTORIES(SYSTEM
                    ${Trilinos_INCLUDE_DIRS}
                    ${Trilinos_TPL_INCLUDE_DIRS})

## Get the absolute path of Trilinos_LIBRARY_DIRS. This is to avoid
## conflicts with Trilinos_TPL_LIBRARY_DIRS which is typically given
## as absolute paths.
## To be removed once Sandia bug #5692, cf.
##   https://software.sandia.gov/bugzilla/show_bug.cgi?id=5692,
## has been resolved.
#GET_FILENAME_COMPONENT(Trilinos_LIBRARY_DIRS_ABSOLUTE
#                       ${Trilinos_LIBRARY_DIRS}
#                       REALPATH)
#
#LINK_DIRECTORIES(${Boost_LIBRARY_DIRS}
#                 ${Trilinos_LIBRARY_DIRS_ABSOLUTE})
# ------------------------------------------------------------------------------
SET(MY_EXECUTABLE "noshExamples")
ADD_EXECUTABLE( ${MY_EXECUTABLE} "exampleCalls.cpp")
TARGET_LINK_LIBRARIES(${MY_EXECUTABLE}
                      "nosh")

# Don't install this example code.
# ------------------------------------------------------------------------------
