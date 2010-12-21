#!/bin/bash

basename="$1"

output="$basename-balanced.nemI"
nem_slice.exe -v -o "$output" -e -m mesh=1x1 -l multikl "$basename.e"

mkdir tmp1
echo "
Input FEM file          = $basename.e
LB file                 = $output
Debug                   = 1
Restart Time list       = off
Reserve space           = nodal=1, elemental=0, global=0
Parallel Disk Info = number=1
Parallel file location = root=tmp,subdir=..
" > nem_spread.inp

nem_spread.exe
rmdir tmp1
