CMAKE_PREFIX_PATH=/opt/trilinos/dev/master/openmpi/1.4.3/gcc/4.7.1/debug/shared:$CMAKE_PREFIX_PATH \
cmake \
    -DCMAKE_INSTALL_PREFIX:PATH=/opt/nosh/dev/master/openmpi/1.4.3/gcc/4.7.2/debug/shared/ \
    -DCMAKE_BUILD_TYPE=Debug \
    -DCMAKE_CXX_COMPILER:STRING=mpic++ \
    -DCMAKE_Fortran_COMPILER:STRING=mpif90 \
    -DBUILD_SHARED_LIBS:BOOL=ON \
    ../../source/
