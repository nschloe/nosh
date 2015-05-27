INCLUDE_DIRECTORIES(${Nosh_SOURCE_DIR})
INCLUDE_DIRECTORIES(
  SYSTEM
  ${Trilinos_INCLUDE_DIRS}
  ${Trilinos_TPL_INCLUDE_DIRS}
  ${EIGEN3_INCLUDE_DIR}
  )

# Get the absolute path of Trilinos_LIBRARY_DIRS. This is to avoid
# conflicts with Trilinos_TPL_LIBRARY_DIRS which is typically given
# as absolute paths.
# To be removed once Sandia bug #5692, cf.
#   https://software.sandia.gov/bugzilla/show_bug.cgi?id=5692,
# has been resolved.
#GET_FILENAME_COMPONENT(Trilinos_LIBRARY_DIRS_ABSOLUTE
#                       ${Trilinos_LIBRARY_DIRS}
#                       REALPATH)

#LINK_DIRECTORIES(${Boost_LIBRARY_DIRS}
#LINK_DIRECTORIES(${Trilinos_LIBRARY_DIRS_ABSOLUTE}
#                 ${Trilinos_TPL_LIBRARY_DIRS})
# ------------------------------------------------------------------------------
SET(MY_EXECUTABLE "nosh-cont")
ADD_EXECUTABLE(${MY_EXECUTABLE} "nosh-cont.cpp" "MyScalarField.cpp")
TARGET_LINK_LIBRARIES(${MY_EXECUTABLE}
                      "nosh")

#SET_TARGET_PROPERTIES(
#  ${MY_EXECUTABLE}
#  PROPERTIES
#  LINKER_LANGUAGE Fortran)

# For now, we don't need COMPONENT installations, cf.
# http://www.cmake.org/pipermail/cmake/2006-October/011362.html.
INSTALL(TARGETS ${MY_EXECUTABLE}
#        # IMPORTANT: Add the executable to the "export-set"
#        EXPORT NoshLibraryDepends
        DESTINATION "${INSTALL_BIN_DIR}")
# ------------------------------------------------------------------------------
