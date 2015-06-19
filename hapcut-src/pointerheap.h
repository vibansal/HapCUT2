
#ifndef _POINTERHEAP_H
#define _POINTERHEAP_H
#include<stdio.h>
#include <math.h>
#include<stdlib.h>
#include<time.h>
#include<string.h>
#include<ctype.h>
#include "common.h"
// binary heap that implements priority queue 
// heap basics: leftchild = 2i+1  rightchild = 2i + 2, parent = i-1/2

// 1 3 6 5 9 8 heap elements  // 0 1 2 3 4 5 array index 

typedef struct PHEAP {
    int length; //int maxlength; 
    int* elements; // where is the element in the array that was originally at index i
} PHEAP;

void pswap(PHEAP* heap, int i, int j);

// trickledown can also be used to update the PHEAP if the score of a node is decreased via an update 
void pmaxHeapify(PHEAP* heap, int node, struct SNPfrags* snpfrag, int* slist);

void pbubbleUp(PHEAP* heap, int node, struct SNPfrags* snpfrag, int* slist); // score of a node is increased, check if the node needs to be bubbled up the heap 
void pbuildmaxheap(PHEAP* heap, struct SNPfrags* snpfrag, int* slist);

int premovemax(PHEAP* heap, struct SNPfrags* snpfrag, int* slist);

void pinitheap(PHEAP* heap, int size);

#endif
