OPT="-O3"
cargo build --release
gcc $OPT -I. -W -Wall -c opusdec.c -o opusdec.o
gcc $OPT -Wall -o opusdec opusdec.o \
	-logg -lm -lspeexdsp target/release/librumpus.a \
	-ldl -lrt -lpthread -lgcc_s -lc -lutil
