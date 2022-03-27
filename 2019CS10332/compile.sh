#!/bin/bash
mpic++ -std=c++11 -Ofast hnsw.cpp -o hnsw -fopenmp
g++ -O2 convert.cpp -o convert
