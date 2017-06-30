all: checkpoint_lib heat_2d

heat_2d: heat-2d.c
	mpicc -g -Wall -o heat_2d heat-2d.c ulcp_files_compress_diff.o ulcp_files_routine.o ulcp_lib.o ulcp_utils.o -lm -std=c99 -lz

checkpoint_lib: 
	mpicc -g -Wall -c ../checkpoint_lib/ulcp_files_compress_diff.c ../checkpoint_lib/ulcp_files_routine.c ../checkpoint_lib/ulcp_lib.c ../checkpoint_lib/ulcp_utils.c

snapshot_integrity: ../checkpoint_lib/snapshot_integrity.cpp ../checkpoint_lib/snapshot_integrity.h ../checkpoint_lib/checkpoint_lib.h
	g++ -g -Wall -o snapshot_integrity ../checkpoint_lib/snapshot_integrity.cpp -std=c++14

clean:
	rm -f heat_2d snapshot_integrity *.o
