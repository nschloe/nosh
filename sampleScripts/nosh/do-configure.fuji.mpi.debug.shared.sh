#CMAKE_PREFIX_PATH=/opt/trilinos/dev/openmpi/1.4.5/gcc/4.7.2/debug/shared:$CMAKE_PREFIX_PATH \
cmake \
    -DCMAKE_INSTALL_PREFIX:PATH=/opt/nosh/dev/openmpi/1.4.5/gcc/4.7.2/debug/shared/ \
    -DCMAKE_BUILD_TYPE=Debug \
    -DCMAKE_CXX_COMPILER:STRING=mpicxx \
    -DCMAKE_Fortran_COMPILER:STRING=mpif90 \
    -DBUILD_SHARED_LIBS:BOOL=ON \
    -DCMAKE_SHARED_LINKER_FLAGS:STRING="-Wl,--no-undefined" \
    ../../source/
