#!/bin/bash
g++ -O3 -c special_c.cpp
gfortran -c special_mod.f90 main.f90 
gfortran main.o special_c.o -lstdc++ -o main
./main

