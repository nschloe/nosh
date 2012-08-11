CMAKE_PREFIX_PATH=/opt/trilinos/dev/master/openmpi/1.4.3/gcc/4.6.3/debug:$CMAKE_PREFIX_PATH \
cmake \
    -DCMAKE_INSTALL_PREFIX:PATH=/opt/nosh/dev/master/openmpi/1.4.3/gcc/4.6.3/debug/ \
    -DCMAKE_BUILD_TYPE=Debug \
    ../../source/
