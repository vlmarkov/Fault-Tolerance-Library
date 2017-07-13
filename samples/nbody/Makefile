all: nbody

nbody: nbody.c
	mpicc -g -Wall -o nbody nbody.c -lm -std=c99

clean:
	rm -f nbody
