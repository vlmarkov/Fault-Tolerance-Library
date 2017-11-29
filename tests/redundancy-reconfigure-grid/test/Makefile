NAME := heat-2d-ulfm
ELF := $(NAME)
SRC := $(NAME).c

CC := mpicc
CFLAGS := -std=c99 -Wall -O2
LDFLAGS := -lm

all: heat-2d

heat-2d: utils grid
	mkdir -p build
	cp slurm.job build/slurm.job
	cp torque.job build/torque.job
	$(CC) -g -Wall -o build/$(ELF) $(SRC) build/utils.o build/grid-task.o $(CFLAGS) $(LDFLAGS)

utils:
	mkdir -p build
	$(CC) -g -Wall -c -o build/utils.o utils.c $(CFLAGS) $(LDFLAGS)

grid:
	mkdir -p build
	$(CC) -g -Wall -c -o build/grid-task.o grid-task.c $(CFLAGS) $(LDFLAGS)

clean:
	rm -Rf build
