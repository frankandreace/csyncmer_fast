CC = g++
CFLAGS = -std=c++17
OPTIMIZE_FLAGS = -march=native -O3 -mavx2 -funroll-loops -fprefetch-loop-arrays
DEBUG_FLAGS = -DDEBUG -g -fno-omit-frame-pointer -Wall -Wextra

# Paths relative to this directory
ROOT_DIR = ../..
SYNG_DIR = ../syng

INCLUDES = -I${ROOT_DIR} -I. -I${SYNG_DIR}
LIBS = -lz -lnthash

# Output
BENCHMARK = benchmark

# Object files
OBJECTS = benchmark.o seqhash.o utils_d.o


### RULES ###

all: CFLAGS += ${OPTIMIZE_FLAGS}
all: clean ${BENCHMARK}

${BENCHMARK}: ${OBJECTS}
	${CC} ${CFLAGS} -o $@ $^ ${LIBS}

benchmark.o: benchmark.cpp
	${CC} ${CFLAGS} ${INCLUDES} -c $< -o $@

seqhash.o: ${SYNG_DIR}/seqhash.c
	${CC} ${CFLAGS} ${INCLUDES} -c $< -o $@

utils_d.o: ${SYNG_DIR}/utils_d.c
	${CC} ${CFLAGS} ${INCLUDES} -c $< -o $@

test: all
	./${BENCHMARK} test

clean:
	rm -f *.o ${BENCHMARK}

debug: CFLAGS += ${DEBUG_FLAGS}
debug: all

.PHONY: all clean debug test
