#!/bin/bash

top_dir=$(pwd)
export RS_SOURCE_DIR=$top_dir

if [ -d "build" ]
then
    echo "deleting existing build directory, starting from scratch (it's ok)"
    rm -R "build"
fi


mkdir -p build &&
    cd build &&

cmake $top_dir &&
    make -B -j4 &&
    make install &&

echo "compiling an example program..."
cd ../example
g++ -std=c++0x scatterFromCascade.cc -o scatterFromCascade `root-config --cflags --glibs --libs` -lRadioScatter
echo "done."
