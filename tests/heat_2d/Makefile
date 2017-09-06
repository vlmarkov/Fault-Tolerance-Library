all: ulcp_lib heat_2d

heat_2d: heat-2d.c
	mpicc -g -Wall -o heat_2d heat-2d.c ulcp.a -L../../3party-libs/zfp-0.5.1/lib/ -lm -lz -lzfp -std=c99 -D ULFM_SUPPORT -D ULFM_TEST

ulcp_lib:
	cd ../../src/ulcp_lib && make all
	mv ../../src/ulcp_lib/ulcp.a ulcp.a

clean:
	rm -f heat_2d *.o *.a
