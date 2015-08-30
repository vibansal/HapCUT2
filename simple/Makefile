CC=gcc -g -D_GNU_SOURCE #-Wall
CFLAGS=-c -Wall

all:SIMPLE

SIMPLE: readinputfiles.o printhaplotypes.c simple.c
	$(CC) readinputfiles.o -o SIMPLE -lm simple.c

readinputfiles.o: readinputfiles.h readinputfiles.c
	$(CC) -c readinputfiles.c

clean:	
	rm -f readinputfiles.o SIMPLE
