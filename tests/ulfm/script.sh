#!/bin/bash

ULFM_PATH=/home/students/mg165/markov/ulfm_src/ulfm_install

export PATH=$ULFM_PATH/bin:$PATH
export LD_LIBRARY_PATH=$ULFM_PATH/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$ULFM_PATH/include:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$ULFM_PATH/glibc:$LD_LIBRARY_PATH

rm -f test

$ULFM_PATH/bin/mpicc -g -Wall -o test test.c -lm -std=c99
