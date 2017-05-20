all: checkpoint_lib heat_2d

heat_2d: heat-2d.c
	mpicc -g -Wall -o heat_2d heat-2d.c checkpoint_lib.o -lm -std=c99 -lz

checkpoint_lib: ../checkpoint_lib/checkpoint_lib.c ../checkpoint_lib/checkpoint_lib.h
	mpicc -c ../checkpoint_lib/checkpoint_lib.c

snapshot_integrity: ../checkpoint_lib/snapshot_integrity.cpp ../checkpoint_lib/snapshot_integrity.h ../checkpoint_lib/checkpoint_lib.h
	g++ -g -Wall -o snapshot_integrity ../checkpoint_lib/snapshot_integrity.cpp -std=c++14

clean:
	rm -f heat_2d snapshot_integrity *.o
