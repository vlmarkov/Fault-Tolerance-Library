NAME := heat-2d-ulfm
ELF := $(NAME)
SRC := $(NAME).c

CC := mpicc
CFLAGS := -std=c99 -Wall -O2
LDFLAGS := -lm

all: heat-2d

heat-2d:
	mkdir -p build
	cp slurm.job build/slurm.job
	cp torque.job build/torque.job
	$(CC) -g -Wall -o build/$(ELF) $(SRC) $(CFLAGS) $(LDFLAGS)

clean:
	rm -Rf build
