CC = g++
CFLAGS = -std=c++17 -march=native -O3 -mavx2 -funroll-loops -fprefetch-loop-arrays
DEBUG_EXTRA = -DDEBUG -g -fno-omit-frame-pointer -Wall -Wextra

# Paths relative to this directory
ROOT_DIR = ../..
SYNG_DIR = ../syng
SIMD_DIR = ../simd

INCLUDES = -I${ROOT_DIR} -I. -I${SYNG_DIR} -I${SIMD_DIR}
LIBS = -lz -lnthash

# Output
BENCHMARK = benchmark
TEST = test

# Shared objects
SHARED_OBJECTS = seqhash.o utils_d.o


### RULES ###

all: clean ${BENCHMARK} ${TEST}

${BENCHMARK}: benchmark.o ${SHARED_OBJECTS}
	${CC} ${CFLAGS} -o $@ $^ ${LIBS}

${TEST}: test.o ${SHARED_OBJECTS}
	${CC} ${CFLAGS} -o $@ $^ ${LIBS}

benchmark.o: benchmark.cpp fasta_reader.h
	${CC} ${CFLAGS} ${INCLUDES} -c $< -o $@

test.o: test.cpp fasta_reader.h
	${CC} ${CFLAGS} ${INCLUDES} -c $< -o $@

seqhash.o: ${SYNG_DIR}/seqhash.c
	${CC} ${CFLAGS} ${INCLUDES} -c $< -o $@

utils_d.o: ${SYNG_DIR}/utils_d.c
	${CC} ${CFLAGS} ${INCLUDES} -c $< -o $@

clean:
	rm -f *.o ${BENCHMARK} ${TEST}

debug: CFLAGS += ${DEBUG_EXTRA}
debug: all

.PHONY: all clean debug
