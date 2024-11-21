#!/bin/bash

if [ ! -d "build" ]; then
    mkdir build
    cmake -S . -B build
fi

cmake --build build || { echo "Build failed"; exit 1; }

time ./build/main || { echo "Program execution failed"; exit 1; }

/usr/local/bin/python3 graph_builder.py || { echo "Graph building script failed"; exit 1; }
