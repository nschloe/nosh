CMAKE_PREFIX_PATH=/opt/trilinos/dev/openmpi/1.4.3/gcc/4.7.1/release/shared:$CMAKE_PREFIX_PATH \
cmake \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX:PATH=/opt/nosh/dev/master/openmpi/1.4.3/gcc/4.7.1/release \
    ../../source/
#    -DBUILD_SHARED_LIBS:BOOL=TRUE \
