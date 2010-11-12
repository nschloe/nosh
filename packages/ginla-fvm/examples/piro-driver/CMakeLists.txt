CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

FIND_PACKAGE( NOX REQUIRED )

FIND_PACKAGE( Boost COMPONENTS system filesystem REQUIRED )
INCLUDE_DIRECTORIES ( ${Boost_INCLUDE_DIRS} )

INCLUDE_DIRECTORIES ( ${ginzburg-landau_SOURCE_DIR}/packages/ginla-fvm/src
                      ${Trilinos_INCLUDE_DIRS}
                      ${Trilinos_TPL_INCLUDE_DIRS} )

LINK_DIRECTORIES( ${Boost_LIBRARY_DIRS}
                  ${Trilinos_LIBRARY_DIRS}
                  ${Trilinos_TPL_LIBRARY_DIRS} )
# ------------------------------------------------------------------------------
SET ( pirodriver_SRCS
      piro-driver.cpp  )
SET ( pirodriver_LIBS
      ginlafdm
      ginla
      komplex2
      recti
      vio
    )
SET( MY_EXECUTABLE
     "piro-driver.exe" )
ADD_EXECUTABLE( ${MY_EXECUTABLE}
                ${pirodriver_SRCS} )
SET_TARGET_PROPERTIES( ${MY_EXECUTABLE}
                       PROPERTIES LINKER_LANGUAGE Fortran )
TARGET_LINK_LIBRARIES( ${MY_EXECUTABLE}
                       ${pirodriver_LIBS}
                       ${Trilinos_LIBRARIES}
                       ${Trilinos_TPL_LIBRARIES}
                       ${Boost_LIBRARIES} )
# ------------------------------------------------------------------------------
SET ( pirodriverfvm_SRCS
      piro-driver-fvm.cpp  )
SET ( pirodriverfvm_LIBS
      ginlafvm-tpetra
      ginla
      komplex2
      vio
    )
SET( MY_EXECUTABLE
     "piro-driver-fvm.exe" )
ADD_EXECUTABLE( ${MY_EXECUTABLE}
                ${pirodriverfvm_SRCS} )
SET_TARGET_PROPERTIES( ${MY_EXECUTABLE}
                       PROPERTIES LINKER_LANGUAGE Fortran )
TARGET_LINK_LIBRARIES( ${MY_EXECUTABLE}
                       ${pirodriverfvm_LIBS}
                       ${Trilinos_LIBRARIES}
                       ${Trilinos_TPL_LIBRARIES}
                       ${Boost_LIBRARIES} )
# ------------------------------------------------------------------------------
SET ( pirodriverfvmepetra_SRCS
      piro-driver-fvm-epetra.cpp  )
SET ( pirodriverfvmepetra_LIBS
      ginlafvm-epetra
      ginla
      komplex2
      vio
    )
SET( MY_EXECUTABLE
     "piro-driver-fvm-epetra.exe" )
ADD_EXECUTABLE( ${MY_EXECUTABLE}
                ${pirodriverfvmepetra_SRCS} )
SET_TARGET_PROPERTIES( ${MY_EXECUTABLE}
                       PROPERTIES LINKER_LANGUAGE Fortran )
TARGET_LINK_LIBRARIES( ${MY_EXECUTABLE}
                       ${pirodriverfvmepetra_LIBS}
                       ${Trilinos_LIBRARIES}
                       ${Trilinos_TPL_LIBRARIES}
                       ${Boost_LIBRARIES} )
# ------------------------------------------------------------------------------
