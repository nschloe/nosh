## http://www.vtk.org/Wiki/VTK/Build_System_Migration
#FIND_PACKAGE(VTK COMPONENTS vtkIOXdmf3 vtkParallelMPI NO_MODULE)
#INCLUDE(${VTK_USE_FILE})

# Treat Trilinos as system code to prevent warnings from Trilinos to appear
# here. Also add Trilinos TPL include directories. Templated TPL packages such
# as Eigen may need this.
INCLUDE_DIRECTORIES(
  SYSTEM
  ${Trilinos_INCLUDE_DIRS}
  ${Trilinos_TPL_INCLUDE_DIRS}
  ${EIGEN3_INCLUDE_DIR}
  #${VTK_INCLUDE_DIRS}
)

FILE(GLOB nosh_SRCS "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")
FILE(GLOB nosh_HEADERS "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp")

ADD_LIBRARY(
  nosh
  ${nosh_SRCS}
  )

# Add target link libraries.
# This has no effect for static builds, but appends the Trilinos libraries to
# the link line to avoid underlinking.  Not adding TARGET_LINK_LIBRARIES()
# makes the link step fail with
# -DCMAKE_SHARED_LINKER_FLAGS:STRING="-Wl,--no-undefined".  We don't need the
# TPL libraries here since the linking information about them is already
# included in the ${Trilinos_LIBRARIES}.  One drawback of this is that the
# development versions of the TPLs must be present on the system since
# Trilinos, by default, doesn't (shared-)link against the fully qualified
# library (e.g., libblas.so.3) but rather its development counterpart (e.g.,
# libblas.so).
TARGET_LINK_LIBRARIES(
  nosh
  ${Trilinos_LIBRARIES}
  ${VTK_LIBRARIES}
  MOAB
  )

SET_TARGET_PROPERTIES(
  nosh PROPERTIES
  PUBLIC_HEADER "${nosh_HEADERS}"
  VERSION ${${PROJECT_NAME}_VERSION}
  SOVERSION ${${PROJECT_NAME}_MAJOR_VERSION}
  )

# Installation of the program
INSTALL(
  TARGETS nosh
  # IMPORTANT: Add an the export (and install it below). It will contain
  # important build information.
  EXPORT NoshLibraryDepends
  # We don't actually need the RUNTIME DESTINATION here since we're only
  # installing a library.
  RUNTIME DESTINATION "${INSTALL_BIN_DIR}" COMPONENT bin
  LIBRARY DESTINATION "${INSTALL_LIB_DIR}" COMPONENT shlib # .so
  ARCHIVE DESTINATION "${INSTALL_LIB_DIR}" COMPONENT shlib # .a
  PUBLIC_HEADER DESTINATION "${INSTALL_INCLUDE_DIR}" COMPONENT dev
  #PUBLIC_HEADER DESTINATION "${INSTALL_INCLUDE_DIR}/nosh" COMPONENT dev
  )

# Install the export set.
INSTALL(
  EXPORT NoshLibraryDepends
  DESTINATION "${INSTALL_CMAKE_DIR}"
  )
