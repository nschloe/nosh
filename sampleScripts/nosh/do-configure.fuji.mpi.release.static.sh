CMAKE_PREFIX_PATH=/opt/trilinos/dev/openmpi/1.4.5/gcc/4.7.2/release/static:$CMAKE_PREFIX_PATH \
cmake \
    -DCMAKE_BUILD_TYPE:STRING=Release \
    -DCMAKE_CXX_COMPILER:STRING=mpic++ \
    -DCMAKE_Fortran_COMPILER:STRING=mpif90 \
    -DCMAKE_INSTALL_PREFIX:PATH=/opt/nosh/dev/master/openmpi/1.4.5/gcc/4.7.2/release/static \
    ../../source/
#    -DBUILD_SHARED_LIBS:BOOL=TRUE \
