#!/bin/bash
set -e   # stop immediately if any command fails

mkdir -p obj mod bin

gfortran -c src/accuracy.f90         -J mod -o obj/accuracy.o
gfortran -c src/error_handling.f90   -J mod -o obj/error_handling.o
gfortran -c src/variables.f90        -J mod -o obj/variables.o
gfortran -c src/setup.f90            -J mod -o obj/setup.o
gfortran -c src/state_vector.f90     -J mod -o obj/state_vector.o
gfortran -c tests/test_state_vector.f90 -J mod -o obj/test_state_vector.o

gfortran obj/accuracy.o obj/error_handling.o obj/variables.o \
         obj/setup.o obj/state_vector.o obj/test_state_vector.o \
         -o bin/test_state_vector

echo "Build OK. Run with: ./bin/test_state_vector"
