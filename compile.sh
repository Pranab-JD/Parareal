#!/bin/bash

### Compile using g++
g++ ../main.cpp -O3 -fopenmp -o Parareal

nvcc ../main.cu -O3 -o Parareal_CUDA -lcublas