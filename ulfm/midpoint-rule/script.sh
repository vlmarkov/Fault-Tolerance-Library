#!/bin/bash

export PATH=/opt/ulfm/bin:$PATH
export LD_LIBRARY_PATH=/opt/ulfm/lib:$LD_LIBRARY_PATH

make clean && make
