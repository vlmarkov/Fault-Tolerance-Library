NAME := heat-2d
ELF := $(NAME)
SRC := $(NAME)

CXX     := mpicxx
CFLAGS  := -g -Wall -O2
LDFLAGS := -lm

all: compute-redundancy

data-redundancy: utils grid
	mkdir -p build
	cp slurm.job build/slurm.job
	cp torque.job build/torque.job
	$(CXX) -o build/$(ELF)-data-redundancy \
	$(SRC)-data-redundancy.cpp build/utils.o \
	build/grid-task.o $(CFLAGS) $(LDFLAGS) -std=c++0x -D MPI_SUPPORT

compute-redundancy: utils grid
	mkdir -p build
	cp slurm.job build/slurm.job
	cp torque.job build/torque.job
	$(CXX) -o build/$(ELF)-compute-redundancy \
	$(SRC)-compute-redundancy.cpp build/utils.o \
	build/grid-task.o $(CFLAGS) $(LDFLAGS) -std=c++0x -D MPI_SUPPORT

utils:
	mkdir -p build
	$(CXX) -c -o build/utils.o utils.cpp \
	$(CFLAGS) $(LDFLAGS) -std=c++0x -D MPI_SUPPORT

grid:
	mkdir -p build
	$(CXX) -c -o build/grid-task.o grid-task.cpp \
	$(CFLAGS) $(LDFLAGS) -std=c++0x -D MPI_SUPPORT

debug:
	mkdir -p build
	g++ $(CFLAGS) -c -o build/debug test.cpp grid-task.cpp \
	$(LDFLAGS) -std=c++0x

clean:
	rm -Rf build
