all: ulcp_lib heat_2d

heat_2d: heat-2d.c
	mpicc -g -Wall -o heat_2d heat-2d.c ulcp.a -lm -lz -std=c99

ulcp_lib:
	cd ../../src/ulcp_lib && make all
	mv ../../src/ulcp_lib/ulcp.a ulcp.a

clean:
	rm -f heat_2d *.o *.a
