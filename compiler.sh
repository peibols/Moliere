#!/bin/bash

g++ -std=c++11 -fpermissive -lgsl -lgslcblas -lm $1.cpp -o $1
