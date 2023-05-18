#!/bin/bash

g++ -std=c++11 -mcmodel=medium -fpermissive -lgsl -lgslcblas -lm $1.cpp -o $1
