OPT="-O0 -fprofile-arcs -ftest-coverage"
cargo build
gcc $OPT -I. -g -W -Wall -c opusdec.c -o opusdec.o
gcc $OPT -g -Wall -o opusdec opusdec.o \
	-logg -lm -lgcov -lspeexdsp target/debug/librumpus.a \
	-ldl -lrt -lpthread -lgcc_s -lc -lutil
