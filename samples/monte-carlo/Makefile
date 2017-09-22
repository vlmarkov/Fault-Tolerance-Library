all: monte-carlo

monte-carlo: monte_carlo.c
	mpicc -g -Wall -o monte_carlo monte_carlo.c -lm -std=c99

clean:
	rm -f monte_carlo
