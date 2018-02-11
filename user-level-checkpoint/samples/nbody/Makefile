all: ulcp_lib nbody

nbody: nbody.c
	mpicc -g -Wall -o nbody nbody.c ulcp.a -lm -lz -std=c99 -D ULFM_SUPPORT

ulcp_lib:
	cd ../../src/ulcp_lib && make all
	mv ../../src/ulcp_lib/ulcp.a ulcp.a

clean:
	rm -f nbody *.o *.a
