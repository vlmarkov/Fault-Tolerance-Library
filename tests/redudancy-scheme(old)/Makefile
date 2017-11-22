all: test

test: grid-task utils
	mkdir -p build
	gcc -g -Wall -o build/test test.c build/grid-task.o build/utils.o -std=c11 -lm

utils:
	mkdir -p build
	gcc -g -Wall -c -o build/utils.o utils.c -std=c11 -lm

grid-task:
	mkdir -p build
	gcc -g -Wall -c -o build/grid-task.o grid-task.c -std=c11 -lm

clean:
	rm -Rf build
