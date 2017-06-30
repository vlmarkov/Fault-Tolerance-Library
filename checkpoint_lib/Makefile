all: ulcp_lib

ulcp_lib: 
	mpicc -g -Wall -c ulcp_files_compress_diff.c ulcp_files_routine.c ulcp_lib.c ulcp_utils.c
	ar cr ulcp.a ulcp_files_compress_diff.o ulcp_files_routine.o ulcp_lib.o ulcp_utils.o

snapshot_integrity: snapshot_integrity.cpp snapshot_integrity.h
	g++ -g -Wall -o snapshot_integrity snapshot_integrity.cpp snapshot_integrity.h -std=c++14

clean:
	rm -f *.a *.o snapshot_integrity
