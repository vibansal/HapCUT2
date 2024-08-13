# HAPCUT2 MAKEFILE
default: all

CFLAGS=-Wall -g -O3 -Wall -D_GNU_SOURCE -I/home/jialechen/software/htslib-1.20/include
LDFLAGS=-L/home/jialechen/software/htslib-1.20/lib -lhts

# DIRECTORIES
B=build
H=hairs-src
X=hapcut2-src
HTSLIB=/home/jialechen/software/htslib-1.20
T=test
# below is the path to CUnit directory, change if need be
CUNIT=/usr/include/CUnit

all: $(B)/extractHAIRS $(B)/HAPCUT2

$(B)/extractHAIRS: $(B)/bamread.o $(B)/hashtable.o $(B)/readvariant.o $(B)/readfasta.o $(B)/hapfragments.o $(H)/extracthairs.c $(H)/parsebamread.c $(H)/realignbamread.c $(H)/nw.c $(H)/realign_pairHMM.c $(H)/estimate_hmm_params.c | $(B)
	$(CC) $(CFLAGS) $(LDFLAGS) -g $(B)/bamread.o $(B)/hapfragments.o $(B)/hashtable.o $(B)/readfasta.o $(B)/readvariant.o -o $(B)/extractHAIRS $(H)/extracthairs.c -pthread -lm -lz -lcurl -llzma -lbz2

$(B)/hapfragments.o: $(H)/hapfragments.c $(H)/hapfragments.h $(H)/readvariant.h | $(B)
	$(CC) -c $(CFLAGS) $(H)/hapfragments.c -o $(B)/hapfragments.o

$(B)/readvariant.o: $(H)/readvariant.c $(H)/readvariant.h $(H)/hashtable.h $(H)/hashtable.c | $(B)
	$(CC) -c $(CFLAGS) $(H)/readvariant.c -o $(B)/readvariant.o

$(B)/bamread.o: $(H)/bamread.h $(H)/bamread.c $(H)/readfasta.h $(H)/readfasta.c | $(B)
	$(CC) -c $(CFLAGS) $(H)/bamread.c -o $(B)/bamread.o

$(B)/hashtable.o: $(H)/hashtable.h $(H)/hashtable.c | $(B)
	$(CC) -c $(CFLAGS) $(H)/hashtable.c -o $(B)/hashtable.o

$(B)/readfasta.o: $(H)/readfasta.c $(H)/readfasta.h | $(B)
	$(CC) -c $(CFLAGS) $(H)/readfasta.c -o $(B)/readfasta.o

# BUILD HAPCUT2

$(B)/HAPCUT2: $(B)/variantgraph.o $(B)/readinputfiles.o $(B)/hapcontig.o $(B)/fragments.o $(B)/readvcf.o $(B)/pointerheap.o $(B)/common.o $(B)/hic.o $(X)/hapcut2.c $(X)/output_phasedvcf.c $(X)/find_maxcut.c $(X)/post_processing.c $(X)/phased_genotyping.c $(X)/maxcut_lr.c | $(B)
	$(CC) $(CFLAGS) $(LDFLAGS) $(B)/common.o $(B)/hic.o $(B)/variantgraph.o $(B)/readinputfiles.o $(B)/hapcontig.o $(B)/fragments.o $(B)/readvcf.o $(B)/pointerheap.o -o $(B)/HAPCUT2 -lm $(X)/hapcut2.c -lhts

$(B)/common.o: $(X)/common.h $(X)/common.c $(X)/variant.h $(X)/fragments.h | $(B)
	$(CC) -c $(CFLAGS) $(X)/common.c -o $(B)/common.o

$(B)/fragments.o: $(X)/fragments.c $(X)/fragments.h $(X)/variant.h $(X)/frag_likelihood.c | $(B)
	$(CC) -c $(CFLAGS) $(X)/fragments.c -o $(B)/fragments.o

$(B)/hic.o: $(X)/hic.h $(X)/hic.c $(X)/common.h | $(B)
	$(CC) -c $(CFLAGS) $(X)/hic.c -o $(B)/hic.o

$(B)/variantgraph.o: $(X)/variantgraph.h $(X)/variantgraph.c $(X)/common.h  | $(B)
	$(CC) -c $(CFLAGS) $(X)/variantgraph.c -o $(B)/variantgraph.o

$(B)/readinputfiles.o: $(X)/readinputfiles.h $(X)/readinputfiles.c $(X)/common.h | $(B)
	$(CC) -c $(CFLAGS) $(X)/readinputfiles.c -o $(B)/readinputfiles.o

$(B)/hapcontig.o: $(X)/hapcontig.c $(X)/common.h | $(B)
	$(CC) -c $(CFLAGS) $(X)/hapcontig.c -o $(B)/hapcontig.o

$(B)/readvcf.o: $(X)/readinputfiles.h $(X)/readvcf.c $(X)/common.h | $(B)
	$(CC) -c $(CFLAGS) $(X)/readvcf.c -o $(B)/readvcf.o

$(B)/pointerheap.o: $(X)/pointerheap.h $(X)/pointerheap.c $(X)/common.h | $(B)
	$(CC) -c $(CFLAGS) $(X)/pointerheap.c -o $(B)/pointerheap.o

# create build directory
$(B):
	mkdir -p $(B)

# INSTALL
install: install-hapcut2 install-hairs

install-hapcut2:
	cp $(B)/HAPCUT2 /usr/local/bin

install-hairs:
	cp $(B)/extractHAIRS /usr/local/bin

# UNINSTALL
uninstall: uninstall-hapcut2 uninstall-hairs

uninstall-hapcut2:
	rm /usr/local/bin/HAPCUT2

uninstall-hairs:
	rm /usr/local/bin/extractHAIRS

# CLEANUP
clean:
	rm -rf $(B)
