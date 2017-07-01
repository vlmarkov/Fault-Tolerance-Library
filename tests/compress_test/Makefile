all: float_point zlib_test

float_point: test_float_point.c
	gcc -g -Wall -o test_float_point test_float_point.c

zlib_test: test_zlib.cpp
	g++ -g -Wall -o test_zlib test_zlib.cpp -std=c++11 -lz

clean:
	rm -f test_float_point test_zlib *.txt
