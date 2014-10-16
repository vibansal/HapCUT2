CC=gcc -Wall -D_GNU_SOURCE
CFLAGS=-c -Wall

HAPCUT: fragmatrix.o readinputfiles.o pointerheap.o hapcut.c find_maxcut.c
	$(CC) fragmatrix.o readinputfiles.o  pointerheap.o -o HAPCUT -lm hapcut.c


HAPCUT-1: fragmatrix.o readinputfiles.o pointerheap.o annealing.o hapcut-annealing.c 
	$(CC) fragmatrix.o readinputfiles.o  pointerheap.o  annealing.o  -o HAPCUT-1 -lm hapcut-annealing.c

all:HAPCUT

fragmatrix.o: fragmatrix.h fragmatrix.c common.h printhaplotypes.c MECscore.c find_starting_haplotypes.c
	$(CC) -c fragmatrix.c

readinputfiles.o: readinputfiles.h readinputfiles.c common.h fragmatrix.h
	$(CC) -c readinputfiles.c

pointerheap.o: pointerheap.h pointerheap.c common.h
	$(CC) -c pointerheap.c

annealing.o: annealing.h annealing.c common.h fragmatrix.h fragmatrix.c
	$(CC) -c annealing.c

clean:	
	rm annealing.o pointerheap.o readinputfiles.o fragmatrix.o HAPCUT
