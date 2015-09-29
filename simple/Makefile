CC=gcc -g -D_GNU_SOURCE -Wall
CFLAGS=-c -Wall

SIMPLE: fragmatrix.o readinputfiles.o simple.c
	$(CC) fragmatrix.o readinputfiles.o -o SIMPLE -lm simple.c

all:SIMPLE

fragmatrix.o: fragmatrix.h fragmatrix.c common.h printhaplotypes.c
	$(CC) -c fragmatrix.c

readinputfiles.o: readinputfiles.h readinputfiles.c common.h fragmatrix.h
	$(CC) -c readinputfiles.c

clean:	
	rm -f readinputfiles.o fragmatrix.o SIMPLE
