
#CC=gcc -Wall
CC=gcc -Wall -D_GNU_SOURCE
CFLAGS=-c -Wall
SAMTOOLS=samtools

all:	hairs


hairs: bamread.o hashtable.o readvariant.o readfasta.o hapfragments.o extracthairs.c 
	$(CC) -I$(SAMTOOLS) -g -O2 bamread.o hapfragments.o hashtable.o readfasta.o readvariant.o -o extractHAIRS extracthairs.c  -L$(SAMTOOLS) -lbam -lm -lz


hapfragments.o:	hapfragments.c hapfragments.h readvariant.h
	$(CC) -c hapfragments.c

readvariant.o: readvariant.c readvariant.h hashtable.h hashtable.c
	$(CC) -c readvariant.c 

bamread.o:	bamread.h bamread.c readfasta.h readfasta.c
	$(CC) -I$(SAMTOOLS) -c bamread.c

hashtable.o: hashtable.h hashtable.c
	$(CC) -c hashtable.c

readfasta.o: readfasta.c readfasta.h
	$(CC) -c readfasta.c

clean:
	rm -f bamread.o readfasta.o readvariant.o hapfragment.o hashtable.o extractHAIRS
