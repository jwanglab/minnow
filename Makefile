CC=gcc
CFLAGS=-std=c99 -O2 # -g for line #s in valgrind

OBJECTS = minnow

all: $(OBJECTS)

minnow: main.c
	$(CC) $(CFLAGS) main.c paf.c incl/minimap2/libminimap2.a -o minnow -lz -lm -lpthread

.PHONY: clean
clean:
	-rm $(OBJECTS)
