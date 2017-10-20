MPICC := mpic++
CFLAGS := -g -Wall

all: heat-2d

heat-2d: utils ft
	$(MPICC) $(CFLAGS) -o heat-2d heat-2d.cpp utils.o fault-tollerance.o

utils:
	$(MPICC) $(CFLAGS) -c utils.cpp

ft:
	$(MPICC) $(CFLAGS) -c fault-tollerance.cpp

clean:
	rm -f heat-2d *.o