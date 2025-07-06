CC = g++ 					# I use this command to specify which compiler I want to use
CFLAGS = -std=c++17	# -march=native -O3 -mavx2 -funroll-loops -fprefetch-loop-arrays # These are 
OPTIMIZE_FLAGS = -march=native -O3 -mavx2 -funroll-loops -fprefetch-loop-arrays
#compilation flags to optimize (max: O3) for the CPU of the machine (native) 
DEBUG_FLAGS = -DDEBUG -g -fno-omit-frame-pointer -Wall -Wextra	# Debug flags

# DIRECTORIES USED
INCLUDE_DIR = src/csyncmer_fast/c_lib/include
DEVELOP_DIR = src/csyncmer_fast/c_lib/dev
EXTERNAL_DIR = src/csyncmer_fast/c_lib/vendor
BUILD_DIR = build
BINARY_DIR = ${BUILD_DIR}/bin
OBJECT_DIR = ${BUILD_DIR}/obj

TESTS_DIR = tests
C_TEST_DIR = ${TESTS_DIR}/c_tests
PY_TEST_DIR = ${TESTS_DIR}/python_tests

# INCLUDING PATHS OF ALL HEADERS FOR COMPILATION
INCLUDES = -I${INCLUDE_DIR} -I${DEVELOP_DIR} -I${EXTERNAL_DIR}

# EXTERNAL LIBRARIES TO BE INCLUDED
C_EXTERNAL_LIBS = -lz -lnthash		# linking with zlib

# EXECUTABLES
C_SPEED_BENCHMARK = ${BINARY_DIR}/speed_benchmark 
C_UNIT_TEST = ${BINARY_DIR}/unit_test

BENCHMARK_SCRIPT = test.sh

# SOURCE FILES FOR THE TEST EXECUTABLE (.c files needed)
C_SPEED_BENCHMARK_SOURCES =  ${C_TEST_DIR}/speed_benchmark.c \
				${EXTERNAL_DIR}/syng/seqhash.c \
				${EXTERNAL_DIR}/syng/utils_d.c

# OBJECT FILES IN THE BUILD DIRECTORY FROM .c FILES
C_SPEED_BENCHMARK_OBJECTS =  ${OBJECT_DIR}/speed_benchmark.o \
				${OBJECT_DIR}/seqhash.o \
				${OBJECT_DIR}/utils_d.o


PYTHON = python3
PYBIND11_INC = $(shell $(PYTHON) -m pybind11 --includes)

PYTHON_SOURCES = python/csyncmer_fast/_bindings.cpp 
PYTHON_MODULE = python/csyncmer_fast/_bindings$(shell $(PYTHON)-config --extension-suffix)


### RULES ###

# DEFAULT TARGET FOR MAKE
all: CFLAGS += ${OPTIMIZE_FLAGS}
all : clean mkdir ${C_SPEED_BENCHMARK} # all depends from mkdir (below) and the binary ${EXE1}

#Build executable from object files
${C_SPEED_BENCHMARK}: ${C_SPEED_BENCHMARK_OBJECTS}
	${CC} ${CFLAGS} -o $@ $^ ${C_EXTERNAL_LIBS}

#Compile test.c into test.o (object file)
${OBJECT_DIR}/speed_benchmark.o: ${C_TEST_DIR}/speed_benchmark.c
	${CC} ${CFLAGS} ${INCLUDES} -c $< -o $@

#Complie seqhash.c into seqhash.o (object file)
${OBJECT_DIR}/seqhash.o: ${EXTERNAL_DIR}/syng/seqhash.c
	${CC} ${CFLAGS} ${INCLUDES} -c $< -o $@

#Compile utils_d.c into utils_d.o (object file)
${OBJECT_DIR}/utils_d.o: ${EXTERNAL_DIR}/syng/utils_d.c
	${CC} ${CFLAGS} ${INCLUDES} -c $< -o $@


#Build unit tests
test: mkdir ${C_UNIT_TEST}

${C_UNIT_TEST}: ${OBJECT_DIR}/unit_test.o
	${CC} ${CFLAGS} -o $@ $^ ${C_EXTERNAL_LIBS}

${OBJECT_DIR}/unit_test.o: ${C_TEST_DIR}/unit_test.c
	${CC} ${CFLAGS} ${INCLUDES} -c $< -o $@

# Run UNIT TESTS
run_test: test python
	./${C_UNIT_TEST}
	pytest -v

# RUN SPEED BENCHMARK
run_speed_benchmark: all
	./${BENCHMARK_SCRIPT}

python:
	python -m build
# python: ${PYTHON_MODULE}
# ${PYTHON_MODULE}: ${PYTHON_SOURCES}
# 	${CC} ${CFLAGS} -I${INCLUDE_DIR} ${PYBIND11_INC} -shared -fPIC $^ -lnthash -o $@ -v

#CREATING DIRECTORIES
mkdir:
	mkdir -p ${BINARY_DIR} ${OBJECT_DIR}

#Clean build artifacts
clean:
	rm -rf ${BUILD_DIR}
	rm -rf ${PYTHON_MODULE}

#Debug build
deubg: CFLAGS += ${DEBUG_FLAGS}
debug: all


.PHONY: all clean mkdir debug test python run_speed_benchmark run_test
