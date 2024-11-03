#!/bin/bash

cmake --build build

./build/main

/usr/local/bin/python3 graph_builder.py
