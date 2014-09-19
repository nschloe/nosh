#CMAKE_PREFIX_PATH=/opt/trilinos/private/:$CMAKE_PREFIX_PATH \
cmake \
    -DCMAKE_INSTALL_PREFIX:PATH=/opt/nosh/dev/openmpi/1.4.5/gcc/4.7.3/debug/shared/ \
    -DCMAKE_BUILD_TYPE=Debug \
    -DCMAKE_CXX_COMPILER:STRING=mpicxx \
    -DBUILD_SHARED_LIBS:BOOL=ON \
    -DCMAKE_SHARED_LINKER_FLAGS:STRING="-Wl,--no-undefined" \
    ../../source/
