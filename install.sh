#!/bin/bash

top_dir=$(pwd)
export RS_SOURCE_DIR=$top_dir

if [ -z "$1" ]
then
    mkdir -p slac_build &&
    cd slac_build 
else
    mkdir -p "$1" &&
    cd "$1" 
fi

cmake $top_dir/slac_rf &&
    make -j4 &&
    make install
