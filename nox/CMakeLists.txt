ENABLE_LANGUAGE( CXX )

INCLUDE_DIRECTORIES ( ${ginzburg-landau_SOURCE_DIR}/libgl
                      /usr/local/trilinos/090603a/include/ )

SET ( ginzburg_landau_SRCS
      gl-nox.cpp )


ADD_EXECUTABLE ( gl-nox.out
                 ${ginzburg_landau_SRCS} )

SET ( internal_LIBS
      gl
    )

# ==============================================================================
# Find Trilinos libraries
FIND_LIBRARY( TEUCHOS teuchos
              HINTS /usr/local/trilinos/current/lib
                    /opt/trilinos-090603a/LINUX_MPI_64BIT/lib )
FIND_LIBRARY( EPETRA epetra
              HINTS /usr/local/trilinos/current/lib
                    /opt/trilinos-090603a/LINUX_MPI_64BIT/lib )
FIND_LIBRARY( EPETRAEXT epetraext
              HINTS /usr/local/trilinos/current/lib
                    /opt/trilinos-090603a/LINUX_MPI_64BIT/lib )
FIND_LIBRARY( NOX nox
              HINTS /usr/local/trilinos/current/lib
                    /opt/trilinos-090603a/LINUX_MPI_64BIT/lib )
FIND_LIBRARY( NOXEPETRA noxepetra
              HINTS /usr/local/trilinos/current/lib
                    /opt/trilinos-090603a/LINUX_MPI_64BIT/lib )
FIND_LIBRARY( AZTECOO aztecoo
              HINTS /usr/local/trilinos/current/lib
                    /opt/trilinos-090603a/LINUX_MPI_64BIT/lib )
FIND_LIBRARY( IFPACK ifpack
              HINTS /usr/local/trilinos/current/lib
                    /opt/trilinos-090603a/LINUX_MPI_64BIT/lib )
FIND_LIBRARY( AMESOS amesos
              HINTS /usr/local/trilinos/current/lib
                    /opt/trilinos-090603a/LINUX_MPI_64BIT/lib )
FIND_LIBRARY( ML ml
              HINTS /usr/local/trilinos/current/lib
                    /opt/trilinos-090603a/LINUX_MPI_64BIT/lib )
# ==============================================================================


# ==============================================================================
# Find other external libs
FIND_LIBRARY( SUPERLU superlu_3.0
              HINTS /usr/local/superlu/3.0/lib
                    /opt/superlu-3.0/lib )
FIND_LIBRARY( BLAS goto-blas-gfortran
              HINTS /usr/local/lib
                    /opt/goto-blas-1.26/lib )
FIND_LIBRARY( LAPACK lapack-gfortran
              HINTS /usr/local/lib
                    /opt/lapack-3.2.1/lib )
FIND_LIBRARY( GFORTRAN gfortran
              HINTS /usr/lib/gcc/x86_64-pc-linux-gnu/4.3.3
                    /opt/gcc-4.3.3/lib64 )
# ==============================================================================


# ==============================================================================
# Set executable linking information.
TARGET_LINK_LIBRARIES ( gl-nox.out
			${internal_LIBS}
                        ${NOX}
                        ${NOXEPETRA}
                        ${AZTECOO}
                        ${IFPACK}
                        ${AMESOS}
                        ${SUPERLU}
                        ${ML}
                        ${EPETRAEXT}
                        ${EPETRA}
                        ${TEUCHOS}
                        ${LAPACK}
                        ${BLAS}
                        ${GFORTRAN}
                      )
# ==============================================================================