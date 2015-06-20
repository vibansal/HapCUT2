
#ifndef INC_hashtable_H
#define INC_hashtable_H

#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>

typedef struct keyvalue {
    char* key;
    int value;
    struct keyvalue* next;
} keyvalue;

keyvalue* keypointer;

typedef struct HASHTABLE {
    int htsize; // prime number that is also the size of HASHTABLE
    int* bucketlengths; // length of each bucket initially 0
    keyvalue** blist; // each bucket is a list of (key,value) pairs, HASHTABLE is an array of buckets
} HASHTABLE;

int hashstring(char* str, int htsize);

void init_hashtable(HASHTABLE* ht);

int insert_keyvalue(HASHTABLE* ht, char* key, int slen, int value);

int getindex(HASHTABLE* ht, char* chrom);

#endif

