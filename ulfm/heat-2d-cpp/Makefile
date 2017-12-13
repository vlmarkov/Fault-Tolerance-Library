NAME := heat-2d
ELF := $(NAME)
SRC := $(NAME)

CXX     := mpicxx
CFLAGS  := -g -Wall -O2
LDFLAGS := -lm

all: data-redundancy compute-redundancy

data-redundancy: utils grid
	mkdir -p build
	cp slurm.job build/slurm.job
	cp torque.job build/torque.job
	$(CXX) -o build/$(ELF)-data-redundancy $(SRC)-data-redundancy.cpp build/utils.o build/grid-task.o $(CFLAGS) $(LDFLAGS)

compute-redundancy: utils grid
	mkdir -p build
	cp slurm.job build/slurm.job
	cp torque.job build/torque.job
	$(CXX) -o build/$(ELF)-compute-redundancy $(SRC)-compute-redundancy.cpp build/utils.o build/grid-task.o $(CFLAGS) $(LDFLAGS)

utils:
	mkdir -p build
	$(CXX) -c -o build/utils.o utils.cpp $(CFLAGS) $(LDFLAGS)

grid:
	mkdir -p build
	$(CXX) -c -o build/grid-task.o grid-task.cpp $(CFLAGS) $(LDFLAGS)

clean:
	rm -Rf build
