ENABLE_LANGUAGE( CXX )

SET ( gl_SRCS 
      glException.cpp
      staggeredGrid.cpp
      ioVirtual.cpp
      ioVtk.cpp
      ioVti.cpp
      ioXdmf.cpp
      ioFactory.cpp
      ginzburgLandau.cpp
      glPrePostOperator.cpp
      glSystem.cpp          )

INCLUDE_DIRECTORIES( /usr/local/trilinos/090603a/include/
                     /opt/trilinos-090603a/LINUX_MPI_64BIT/include/
                     /opt/boost/1.39.0/include/boost-1_39/ 
                     /home/nschloe/svn/paracont/trunk/include/ 
                     /opt/hdf5-1.6.9/include
                     /opt/openmpi-1.3.3/include )

ADD_LIBRARY ( gl STATIC
              ${gl_SRCS} )

TARGET_LINK_LIBRARIES ( gl )
