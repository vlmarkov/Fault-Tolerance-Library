NAME := heat-2d
ELF := $(NAME)
SRC := $(NAME)

CXX     := mpicxx
CFLAGS  := -g -Wall -O2
LDFLAGS := -lm

all: compute-redundancy

compute-redundancy: utils grid redundancy
	mkdir -p build
	cp slurm.job build/slurm.job
	cp torque.job build/torque.job
	$(CXX) -o build/$(ELF)-compute-redundancy \
	$(SRC)-compute-redundancy.cpp build/utils.o build/redundancy.o \
	build/grid-task.o $(CFLAGS) $(LDFLAGS) -std=c++0x -D MPI_SUPPORT

utils:
	mkdir -p build
	$(CXX) -c -o build/utils.o src/utils.cpp \
	$(CFLAGS) $(LDFLAGS) -std=c++0x -D MPI_SUPPORT

grid:
	mkdir -p build
	$(CXX) -c -o build/grid-task.o src/grid-task.cpp \
	$(CFLAGS) $(LDFLAGS) -std=c++0x -D MPI_SUPPORT

redundancy:
	mkdir -p build
	$(CXX) -c -o build/redundancy.o src/redundancy.cpp \
	$(CFLAGS) $(LDFLAGS) -std=c++0x -D MPI_SUPPORT

debug:
	mkdir -p build
	g++ $(CFLAGS) -o build/debug test.cpp \
	src/grid-task.cpp \
	src/utils.cpp \
	src/redundancy.cpp \
	$(LDFLAGS) -std=c++0x

clean:
	rm -Rf build
