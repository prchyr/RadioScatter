#!/bin/bash

top_dir=$(pwd)
export RS_SOURCE_DIR=$top_dir

if [ -d "build" ]
then
    echo "deleting existing build directory, starting from scratch (it's ok)"
    rm -R "build"
fi

# if [ -z "$1" ]
# then
#     mkdir -p slac_build &&
#     cd slac_build 
# else
#     mkdir -p "$1" &&
#     cd "$1" 
# fi

mkdir -p build &&
    cd build

cmake $top_dir &&
    make -B -j4 &&
    make install
