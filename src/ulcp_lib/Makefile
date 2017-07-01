all: ulcp_lib

ulcp_lib: ulcp.h utils files_routine files_compress_diff main
	ar cr ulcp.a files_compress_diff.o files_routine.o main.o utils.o
	rm -f *.o

utils: utils.c utils.h
	mpicc -g -Wall -Werror -c utils.c

files_routine: files_routine.c files_routine.h
	mpicc -g -Wall -Werror -c files_routine.c

files_compress_diff: files_compress_diff.c files_compress_diff.h
	mpicc -g -Wall -Werror -c files_compress_diff.c -lz

main: main.c main.h
	mpicc -g -Wall -Werror -c main.c

clean:
	rm -f *.a
