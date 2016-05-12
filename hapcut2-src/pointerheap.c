// binary heap that implements priority queue 
// heap basics: leftchild = 2i+1  rightchild = 2i + 2, parent = i-1/2

// 1 3 6 5 9 8 heap elements
// 0 1 2 3 4 5 array index 
#include "pointerheap.h"

void pswap(PHEAP* heap, int i, int j) {
    //	heap->location[heap->harray[i].index] = i;  heap->location[heap->harray[j].index] = j; 
    // need to maintain location of each SNP in heap 
}

// trickledown can also be used to update the PHEAP if the score of a node is decreased via an update 

void pmaxHeapify(PHEAP* heap, int node, struct SNPfrags* snpfrag, int* slist) {
    int lc = 2 * node + 1;
    int rc = 2 * node + 2;
    int maxindex = node;
    int temp;
    //	printf("inside maxheapify %d %d %d\n",lc,rc,node);
    if (rc >= heap->length) {
        if (lc < heap->length) maxindex = lc;
    } else {
        if (fabsf(snpfrag[slist[heap->elements[lc]]].score) >= fabsf(snpfrag[slist[heap->elements[rc]]].score)) maxindex = lc;
        else maxindex = rc;
    }
    if (fabsf(snpfrag[slist[heap->elements[node]]].score) < fabsf(snpfrag[slist[heap->elements[maxindex]]].score)) // swap and percolate down
    {
        //temp = snpfrag[slist[node]].heaploc; snpfrag[slist[node]].heaploc = snpfrag[slist[maxindex]].heaploc; snpfrag[slist[maxindex]].heaploc = temp;
        temp = heap->elements[node];
        heap->elements[node] = heap->elements[maxindex];
        heap->elements[maxindex] = temp;
        snpfrag[slist[heap->elements[node]]].heaploc = node;
        snpfrag[slist[heap->elements[maxindex]]].heaploc = maxindex;
        pmaxHeapify(heap, maxindex, snpfrag, slist);
    }

}

void pbubbleUp(PHEAP* heap, int node, struct SNPfrags* snpfrag, int* slist) // score of a node is increased, check if the node needs to be bubbled up the heap 
{
    int parent, temp;
    
    while (node > 0) {
        parent = (node - 1);
        parent /= 2;/*
        fprintf(stderr,"spot1\n");
        fprintf(stderr,"%d\n",node);
        fprintf(stderr,"%d\n",heap->elements[node]);
        fprintf(stderr,"%d\n",slist[heap->elements[node]]);
        fprintf(stderr,"spot2\n");*/
        //fprintf(stderr,"%f\n",fabsf(snpfrag[slist[heap->elements[parent]]].score));
        //fprintf(stderr,"spot3");
        //printf("node %d score %f parent %d %f \n",node,heap->harray[node].score,parent,heap->harray[parent].score);
        if (fabsf(snpfrag[slist[heap->elements[node]]].score) > fabsf(snpfrag[slist[heap->elements[parent]]].score)) {
            //temp = snpfrag[slist[node]].heaploc; snpfrag[slist[node]].heaploc = snpfrag[slist[parent]].heaploc; snpfrag[slist[parent]].heaploc = temp;
            temp = heap->elements[node];
            heap->elements[node] = heap->elements[parent];
            heap->elements[parent] = temp;
            snpfrag[slist[heap->elements[node]]].heaploc = node;
            snpfrag[slist[heap->elements[parent]]].heaploc = parent;
        } else break;
        node = parent;
    }
}

void pbuildmaxheap(PHEAP* heap, struct SNPfrags* snpfrag, int* slist) {
    int i = 0;
    for (i = heap->length / 2 - 1; i >= 0; i--) pmaxHeapify(heap, i, snpfrag, slist);
    //fprintf(stdout,"heapify %d hl %d\n",i,heap->length);
}

int premovemax(PHEAP* heap, struct SNPfrags* snpfrag, int* slist) {
    // easy to copy elements of the harray[0] element from the heap 
    if (heap->length == 0) return -1;
    heap->elements[0] = heap->elements[heap->length - 1];
    snpfrag[slist[heap->elements[0]]].heaploc = 0;
    heap->length--;
    if (heap->length == 0) return 1;
    // swap root's value with the maximum of it's two children, maxHeapify 
    pmaxHeapify(heap, 0, snpfrag, slist);
    return 1;

}

void pinitheap(PHEAP* heap, int size) {
    heap->elements = (int*) calloc(size,sizeof (int));
    heap->length = size; //heap->maxlength = size; 
}

/*
int mainchecl (int argc, char** argv)
{
//	float f= -3.4;  float f1 = fabsf(f); printf("abs %f %f \n",f,f1);
        int i =0;
        PHEAP heap;
        createheap(&heap,7);
        heap.harray[0].score = 3; heap.harray[1].score = -50; heap.harray[2].score = 6; heap.harray[3].score = 8; heap.harray[4].score = 9; 
        heap.harray[5].score = 12; heap.harray[6].score = 14;
        for (i=0;i<heap.length;i++) fprintf(stdout,"elem %d %f \n",i,heap.harray[i].score);
        buildmaxheap(&heap);
        for (i=0;i<heap.length;i++) fprintf(stdout,"elem %d %f \n",i,heap.harray[i].score);
        heap.harray[5].score = -553; 
        bubbleUp(&heap,5);
//	maxHeapify(&heap,2);
        for (i=0;i<heap.length;i++) fprintf(stdout,"elem %d %f \n",i,heap.harray[i].score);


        return 1;
}*/
