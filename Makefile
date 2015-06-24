
# HAPCUT MAKEFILE
default: all

CC=gcc -Wall -D_GNU_SOURCE
CFLAGS=-c -Wall

# DIRECTORIES
B=build
H=hairs-src
X=hapcut-src
HTSLIB=submodules/htslib
SAMTOOLS=submodules/samtools
T=test
# below is the path to CUnit directory, change if need be
CUNIT=/usr/include/CUnit


all: $(B)/extractHAIRS $(B)/extractFOSMID $(B)/HAPCUT

# TEST
# requires Cunit
tests:
	$(CC) $(T)/test.c -o $(B)/test -lcunit -I $(CUNIT)
	./$(B)/test

# check if git present, else print error message
# credit to lhunath for the one-line check below http://stackoverflow.com/questions/592620/check-if-a-program-exists-from-a-bash-script
git-exists:
	@hash git 2>/dev/null || { echo >&2 "Git not installed (required for auto-download of required submodules). Install git and retry, or manually download samtools 1.2 and htslib 1.2.1 and place the unzipped folders in \"submodules\" directory with names matching their Makefile variables (default \"samtools\" and \"htslib\")"; exit 1; }

# if samtools makefile not present, then submodules have not yet been downloaded (init & updated)
$(SAMTOOLS)/Makefile: git-exists
	git submodule init
	git submodule update

# just call the samtools makefile target because it downloads both htslib and samtools submodules
$(HTSLIB)/Makefile: $(SAMTOOLS)/Makefile

$(SAMTOOLS)/libbam.a: $(SAMTOOLS)/Makefile
	echo "Building Samtools bam library..."
	make -C $(SAMTOOLS) lib

$(HTSLIB)/libhts.a: $(HTSLIB)/Makefile
	echo "Building HTSlib libraries..."
	make -C $(HTSLIB) lib-static

# BUILD HAIRS

$(B)/extractHAIRS: $(B)/bamread.o $(B)/hashtable.o $(B)/readvariant.o $(B)/readfasta.o $(B)/hapfragments.o $(H)/extracthairs.c $(SAMTOOLS)/libbam.a $(HTSLIB)/libhts.a | $(B)
	$(CC) -I$(SAMTOOLS) -I$(HTSLIB) -g -O2 $(B)/bamread.o $(B)/hapfragments.o $(B)/hashtable.o $(B)/readfasta.o $(B)/readvariant.o -o $(B)/extractHAIRS $(H)/extracthairs.c  -L$(SAMTOOLS) -L$(HTSLIB) -pthread -lhts -lbam -lm -lz

$(B)/extractFOSMID: $(B)/bamread.o $(B)/hashtable.o $(B)/readvariant.o $(B)/readfasta.o $(B)/hapfragments.o $(H)/extracthairs.c $(H)/fosmidbam_hairs.c $(H)/print_clusters.c $(SAMTOOLS)/libbam.a $(HTSLIB)/libhts.a | $(B)
	$(CC) -I$(SAMTOOLS) -I$(HTSLIB) -g -O2 $(B)/bamread.o $(B)/hapfragments.o $(B)/hashtable.o $(B)/readfasta.o $(B)/readvariant.o -o $(B)/extractFOSMID $(H)/extracthairs.c  -L$(SAMTOOLS) -L$(HTSLIB) -pthread -lhts -lbam -lm -lz

#INDELCOUNTS: $(B)/bamread.o $(B)/hashtable.o $(B)/readvariant.o $(B)/readfasta.o $(B)/hapfragments.o $(H)/indelcounts.c | $(B)
#	$(CC) -I$(SAMTOOLS) -g -O2 $(B)/bamread.o $(B)/hapfragments.o $(B)/hashtable.o $(B)/readfasta.o $(B)/readvariant.o -o $(B)/INDELCOUNTS $(H)/indelcounts.c -L$(SAMTOOLS) -lbam -lm -lz

#(B)/calculateGLL: $(B)/bamread.o $(B)/hashtable.o $(B)/readvariant.o $(B)/readfasta.o $(B)/hapfragments.o $(H)/calculateGLL.c | $(B)
#	$(CC) -I$(SAMTOOLS) -g -O2 $(B)/bamread.o $(B)/hapfragments.o $(B)/hashtable.o $(B)/readfasta.o $(B)/readvariant.o -o $(B)/calculateGLL $(H)/calculateGLL.c  -L$(SAMTOOLS) -lbam -lm -lz

#(B)/PCR: $(B)/bamread.o $(B)/hashtable.o $(B)/readvariant.o $(B)/readfasta.o $(B)/hapfragments.o $(H)/calculateGLL.c | $(B)
#	$(CC) -I$(SAMTOOLS) -g -O2 $(B)/bamread.o $(B)/hapfragments.o $(B)/hashtable.o $(B)/readfasta.o $(B)/readvariant.o -o $(B)/PCR $(H)/PCRdups.c  -L$(SAMTOOLS) -lbam -lm -lz

$(B)/hapfragments.o: $(H)/hapfragments.c $(H)/hapfragments.h $(H)/readvariant.h | $(B)
	$(CC) -c $(H)/hapfragments.c -o $(B)/hapfragments.o

$(B)/readvariant.o: $(H)/readvariant.c $(H)/readvariant.h $(H)/hashtable.h $(H)/hashtable.c | $(B)
	$(CC) -c $(H)/readvariant.c -o $(B)/readvariant.o

$(B)/bamread.o: $(H)/bamread.h $(H)/bamread.c $(H)/readfasta.h $(H)/readfasta.c $(SAMTOOLS)/libbam.a $(HTSLIB)/libhts.a | $(B)
	$(CC) -I$(SAMTOOLS) -I$(HTSLIB) -c $(H)/bamread.c -o $(B)/bamread.o

$(B)/hashtable.o: $(H)/hashtable.h $(H)/hashtable.c | $(B)
	$(CC) -c $(H)/hashtable.c -o $(B)/hashtable.o

$(B)/readfasta.o: $(H)/readfasta.c $(H)/readfasta.h | $(B) 
	$(CC) -c $(H)/readfasta.c -o $(B)/readfasta.o

# BUILD HAPCUT

$(B)/HAPCUT: $(B)/fragmatrix.o $(B)/readinputfiles.o $(B)/pointerheap.o $(X)/hapcut.c $(X)/find_maxcut.c | $(B) 
	$(CC) $(B)/fragmatrix.o $(B)/readinputfiles.o $(B)/pointerheap.o -o $(B)/HAPCUT -lm $(X)/hapcut.c

#$(B)/HAPCUT-1: fragmatrix.o readinputfiles.o pointerheap.o annealing.o $(X)/hapcut-annealing.c | $(B)
#	$(CC) $(X)/fragmatrix.o $(X)/readinputfiles.o $(X)/pointerheap.o $(X)/annealing.o  -o $(B)/HAPCUT-1 -lm $(X)/hapcut-annealing.c

$(B)/fragmatrix.o: $(X)/fragmatrix.h $(X)/fragmatrix.c $(X)/common.h $(X)/printhaplotypes.c $(X)/MECscore.c $(X)/find_starting_haplotypes.c | $(B)
	$(CC) -c $(X)/fragmatrix.c -o $(B)/fragmatrix.o

$(B)/readinputfiles.o: $(X)/readinputfiles.h $(X)/readinputfiles.c $(X)/common.h $(X)/fragmatrix.h | $(B)
	$(CC) -c $(X)/readinputfiles.c -o $(B)/readinputfiles.o

$(B)/pointerheap.o: $(X)/pointerheap.h $(X)/pointerheap.c $(X)/common.h | $(B)
	$(CC) -c $(X)/pointerheap.c -o $(B)/pointerheap.o

$(B)/annealing.o: $(X)/annealing.h $(X)/annealing.c $(X)/common.h $(X)/fragmatrix.h $(X)/fragmatrix.c | $(B)
	$(CC) -c $(X)/annealing.c -o $(B)/annealing.o

# create build directory
$(B):
	mkdir -p $(B)

# INSTALL
install: install-hapcut install-hairs install-fosmid

install-hapcut:
	cp $(B)/HAPCUT /bin

install-hairs:
	cp $(B)/extractHAIRS /bin

install-fosmid:
	cp $(B)/extractFOSMID /bin

# UNINSTALL
uninstall: uninstall-hapcut uninstall-hairs uninstall-fosmid

uninstall-hapcut:
	rm /bin/HAPCUT

uninstall-hairs:
	rm /bin/extractHAIRS

uninstall-fosmid:
	rm /bin/extractFOSMID

# CLEANUP
nuke: clean
	make clean -C submodules/samtools
	make clean -C submodules/htslib

clean:
	rm -rf $(B)