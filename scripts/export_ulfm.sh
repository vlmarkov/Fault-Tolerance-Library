#!/bin/bash

export PATH=$HOME/opt/ulfm/bin:$PATH
export LD_LIBRARY_PATH=$HOME/opt/ulfm/lib:$LD_LIBRARY_PATH



shell$ ./configure --prefix=/home/students/mg165/markov/ulfm_target \
       --enable-mpi-ext=ftmpi --with-ft=mpi \
       --disable-io-romio --enable-contrib-no-build=vt \
       --with-platform=optimized \
       CC=gcc CXX=g++ F77=gfortran FC=gfortran


--prefix=/home/students/mg165/markov/ulfm_target


export PATH=/home/students/mg165/markov/ulfm_target/bin:$PATH
export LD_LIBRARY_PATH=/home/students/mg165/markov/ulfm_target/lib:$LD_LIBRARY_PATH

LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/students/mg165/markov/ulfm_target/lib/
export LD_LIBRARY_PATH

./configure --prefix=/home/students/mg165/markov/ulfm_src/ulfm_2_0 \
--with-ft=mpi \
--with-tm=/opt/torque-6.1.1.1