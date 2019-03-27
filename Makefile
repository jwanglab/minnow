CC=gcc
CFLAGS=-std=c99 -O2

OBJECTS = minnow

all: $(OBJECTS)

minnow: main.c
	$(CC) $(CFLAGS) main.c incl/minimap2/libminimap2.a -o minnow -lz -lm -lpthread

.PHONY: clean
clean:
	-rm $(OBJECTS)
