all: checkpoint_lib heat_2d

heat_2d: heat-2d.c
	mpicc -g -Wall -o heat_2d heat-2d.c checkpoint_lib.o -lm -std=c99

checkpoint_lib: ../checkpoint_lib/checkpoint_lib.c ../checkpoint_lib/checkpoint_lib.h
	mpicc -c ../checkpoint_lib/checkpoint_lib.c

clean:
	rm -f heat_2d *.o
