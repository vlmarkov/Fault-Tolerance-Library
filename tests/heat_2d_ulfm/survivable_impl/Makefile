all: heat-2d

heat-2d: utils
	mpicc -g -Wall -O3 -o heat-2d heat-2d.c utils.o -lm

utils: utils.h utils.c
	mpicc -g -Wall -O3 -c utils.c

clean:
	rm -f heat-2d *.o
