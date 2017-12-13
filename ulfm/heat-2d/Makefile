NAME := heat-2d
ELF := $(NAME)
SRC := $(NAME)

CC := mpicc
CFLAGS := -std=c99 -Wall -O2
LDFLAGS := -lm

all: data-redundancy compute-redundancy

data-redundancy: utils grid
	mkdir -p build
	cp slurm.job build/slurm.job
	cp torque.job build/torque.job
	$(CC) -g -Wall -o build/$(ELF)-data-redundancy $(SRC)-data-redundancy.c build/utils.o build/grid-task.o $(CFLAGS) $(LDFLAGS)

compute-redundancy: utils grid
	mkdir -p build
	cp slurm.job build/slurm.job
	cp torque.job build/torque.job
	$(CC) -g -Wall -o build/$(ELF)-compute-redundancy $(SRC)-compute-redundancy.c build/utils.o build/grid-task.o $(CFLAGS) $(LDFLAGS)

utils:
	mkdir -p build
	$(CC) -g -Wall -c -o build/utils.o utils.c $(CFLAGS) $(LDFLAGS)

grid:
	mkdir -p build
	$(CC) -g -Wall -c -o build/grid-task.o grid-task.c $(CFLAGS) $(LDFLAGS)

clean:
	rm -Rf build
