git clone https://github.com/vibansal/HapCUT2.git
cd hapcut2
git clone https://github.com/samtools/htslib.git
cd htslib
autoreconf -i
./configure
make
cd ..
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$(pwd)/htslib"
make
