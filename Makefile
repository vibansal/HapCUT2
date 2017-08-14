
# HAPCUT2 MAKEFILE
default: all

CC=gcc -g -O3 -Wall -D_GNU_SOURCE
CFLAGS=-c -Wall

# DIRECTORIES
B=build
H=hairs-src
X=hapcut2-src
HTSLIB=submodules/htslib
SAMTOOLS=submodules/samtools
T=test
# below is the path to CUnit directory, change if need be
CUNIT=/usr/include/CUnit

all: $(B)/extractHAIRS $(B)/extractFOSMID $(B)/HAPCUT2

# if samtools makefile not present, then submodules have not yet been downloaded (init & updated)
# first check if git present, else print error message
# credit to lhunath for the one-line check for git below http://stackoverflow.com/questions/592620/check-if-a-program-exists-from-a-bash-script
$(SAMTOOLS)/Makefile:
	@hash git 2>/dev/null || { echo >&2 "Git not installed (required for auto-download of required submodules). Install git and retry, or manually download samtools 1.2 and htslib 1.2.1 and place the unzipped folders in \"submodules\" directory with names matching their Makefile variables (default \"samtools\" and \"htslib\")"; exit 1; }
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

#temporarily removed -O2 flag after -I$(HTSLIB)
$(B)/extractHAIRS: $(B)/bamread.o $(B)/hashtable.o $(B)/readvariant.o $(B)/readfasta.o $(B)/hapfragments.o $(H)/extracthairs.c $(SAMTOOLS)/libbam.a $(HTSLIB)/libhts.a $(H)/parsebamread.c $(H)/realignbamread.c $(H)/nw.c | $(B)
	$(CC) -I$(SAMTOOLS) -I$(HTSLIB) -g $(B)/bamread.o $(B)/hapfragments.o $(B)/hashtable.o $(B)/readfasta.o $(B)/readvariant.o -o $(B)/extractHAIRS $(H)/extracthairs.c  -L$(SAMTOOLS) -L$(HTSLIB) -pthread -lhts -lbam -lm -lz
#temporarily removed -O2 flag after -I$(HTSLIB)
$(B)/extractFOSMID: $(B)/bamread.o $(B)/hashtable.o $(B)/readvariant.o $(B)/readfasta.o $(B)/hapfragments.o $(H)/extracthairs.c $(H)/fosmidbam_hairs.c $(H)/print_clusters.c $(SAMTOOLS)/libbam.a $(HTSLIB)/libhts.a | $(B)
	$(CC) -I$(SAMTOOLS) -I$(HTSLIB) -g $(B)/bamread.o $(B)/hapfragments.o $(B)/hashtable.o $(B)/readfasta.o $(B)/readvariant.o -o $(B)/extractFOSMID $(H)/extracthairs.c  -L$(SAMTOOLS) -L$(HTSLIB) -pthread -lhts -lbam -lm -lz

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

# BUILD HAPCUT2

$(B)/HAPCUT2: $(B)/fragmatrix.o $(B)/readinputfiles.o $(B)/pointerheap.o $(B)/common.o $(X)/hapcut2.c $(X)/find_maxcut.c $(X)/post_processing.c| $(B)
	$(CC) $(B)/common.o $(B)/fragmatrix.o $(B)/readinputfiles.o $(B)/pointerheap.o -o $(B)/HAPCUT2 -lm $(X)/hapcut2.c

$(B)/common.o: $(X)/common.h $(X)/common.c | $(B)
	$(CC) -c $(X)/common.c -o $(B)/common.o

$(B)/fragmatrix.o: $(X)/fragmatrix.h $(X)/fragmatrix.c $(X)/common.h $(X)/printhaplotypes.c | $(B)
	$(CC) -c $(X)/fragmatrix.c -o $(B)/fragmatrix.o

$(B)/readinputfiles.o: $(X)/readinputfiles.h $(X)/readinputfiles.c $(X)/common.h $(X)/fragmatrix.h | $(B)
	$(CC) -c $(X)/readinputfiles.c -o $(B)/readinputfiles.o

$(B)/pointerheap.o: $(X)/pointerheap.h $(X)/pointerheap.c $(X)/common.h | $(B)
	$(CC) -c $(X)/pointerheap.c -o $(B)/pointerheap.o

# create build directory
$(B):
	mkdir -p $(B)

# INSTALL
install: install-hapcut2 install-hairs install-fosmid

install-samtools:
	make -C $(SAMTOOLS)
	make -C $(SAMTOOLS) install

install-hapcut2:
	cp $(B)/HAPCUT2 /usr/local/bin

install-hairs:
	cp $(B)/extractHAIRS /usr/local/bin

install-fosmid:
	cp $(B)/extractFOSMID /usr/local/bin

# UNINSTALL
uninstall: uninstall-hapcut2 uninstall-hairs uninstall-fosmid

uninstall-hapcut2:
	rm /usr/local/bin/HAPCUT2

uninstall-hairs:
	rm /usr/local/bin/extractHAIRS

uninstall-fosmid:
	rm /usr/local/bin/extractFOSMID

# CLEANUP
nuke: clean
	make clean -C submodules/samtools
	make clean -C submodules/htslib

clean:
	rm -rf $(B)
