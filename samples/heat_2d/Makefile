heat := heat-2d
heat_objs := heat-2d.o

CC := mpicc
CFLAGS := -std=c99 -Wall -O2
LDFLAGS := -lm

.PHONY: all clean

all: $(heat)

$(heat): $(heat_objs)
	$(CC) -o $@ $^ $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

heat-2d.o: heat-2d.c

clean:
	@rm -rf *.o $(heat)
