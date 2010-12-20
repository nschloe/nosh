#!/bin/bash

nem_slice -e -m mesh=2x2 -l multikl "$1"

#mkdir tmp1
echo "
Input FEM file          = $1
LB file                 = $1-m4-bKL.nemI 
Debug                   = 1
Restart Time list       = off
Reserve space           = nodal=1, elemental=0, global=m
marallel Disk Info = number=1
Parallel file location = root=tmp,subdir=..
" > nem_spread.inp
nem_spread
#rmdir tmp1
