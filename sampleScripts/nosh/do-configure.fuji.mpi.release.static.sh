CMAKE_PREFIX_PATH=/opt/trilinos/dev/openmpi/1.4.3/gcc/4.7.1/release/static:$CMAKE_PREFIX_PATH \
cmake \
    -DCMAKE_BUILD_TYPE:STRING=Release \
    -DCMAKE_CXX_COMPILER:STRING=mpicxx \
    -DCMAKE_INSTALL_PREFIX:PATH=/opt/nosh/dev/master/openmpi/1.4.3/gcc/4.7.1/release/static \
    ../../source/
#    -DBUILD_SHARED_LIBS:BOOL=TRUE \
