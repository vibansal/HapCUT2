#include "hashtable.h"

// this is a basic hashtable that was designed for mapping chromosomes to integers (chrx - 0) 
//  used in extract_hairs for haplotype assembly as well as indel realignment 

// simple hash function that takes a string and returns an integer hash value

int hashstring(char* str, int htsize) {
    unsigned long hash = 5381;
    int c;
    while ((c = *str++)) {
        hash = ((hash << 5) + hash) + c;
        if (hash >= htsize) hash = hash % htsize;
    }
    return hash;
}

void init_hashtable(HASHTABLE* ht) {
    int i = 0;
    ht->bucketlengths = (int*) malloc(sizeof (int)*ht->htsize);
    for (i = 0; i < ht->htsize; i++) ht->bucketlengths[i] = 0;
    ht->blist = (keyvalue**) malloc(sizeof (keyvalue*) * ht->htsize);
    for (i = 0; i < ht->htsize; i++) ht->blist[i] = NULL; //(keyvalue*)malloc(sizeof(keyvalue)*20);
}

int insert_keyvalue(HASHTABLE* ht, char* key, int slen, int value) {
    int hash = hashstring(key, ht->htsize);
    keyvalue* tempkey = (keyvalue*) malloc(sizeof (keyvalue));
    tempkey->value = value;
    tempkey->key = (char*) malloc(slen + 1);
    int i = 0;
    for (i = 0; i < slen; i++) tempkey->key[i] = key[i];
    tempkey->key[i] = '\0';
    tempkey->next = ht->blist[hash];
    ht->blist[hash] = tempkey;
    ht->bucketlengths[hash]++;
    return 1;
}

int getindex(HASHTABLE* ht, char* chrom) {
    int hash = hashstring(chrom, ht->htsize);
    keypointer = ht->blist[hash];
    while (keypointer != NULL) {
        if (strcmp(keypointer->key, chrom) == 0) return keypointer->value;
        keypointer = keypointer->next;
    }
    return -1;
}


