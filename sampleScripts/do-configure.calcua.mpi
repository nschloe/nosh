#!/bin/bash

#module purge
module load Boost/1.43.0-ictce-3.2.1.015.u1
module load VTK/5.6.1-ictce-3.2.1.015.u1

cmake \
    -D Trilinos_DIR:PATH="$VSC_DATA/software/trilinos/dev/master/install/include/" \
    -D BOOST_ROOT:PATH="${SOFTROOTBOOST}" \
    ../../source
