CC = g++ 					# I use this command to specify which compiler I want to use
CFLAGS = -std=c++17	# -march=native -O3 -mavx2 -funroll-loops -fprefetch-loop-arrays # These are 
OPTIMIZE_FLAGS = -march=native -O3 -mavx2 -funroll-loops -fprefetch-loop-arrays
#compilation flags to optimize (max: O3) for the CPU of the machine (native) 
DEBUG_FLAGS = -DDEBUG -g -fno-omit-frame-pointer -Wall -Wextra	# Debug flags

# DIRECTORIES USED
INCLUDE_DIR = include
DEVELOP_DIR = dev
EXTERNAL_DIR = external
TEST_DIR = tests
BUILD_DIR = build
BINARY_DIR = ${BUILD_DIR}/bin
OBJECT_DIR = ${BUILD_DIR}/obj

# INCLUDING PATHS OF ALL HEADERS FOR COMPILATION
INCLUDES = -I${INCLUDE_DIR} -I${DEVELOP_DIR} -I${EXTERNAL_DIR}

# EXTERNAL LIBRARIES TO BE INCLUDED
EXE1_LIBS = -lz -lnthash		# linking with zlib

# EXECUTABLE
EXE1 = ${BINARY_DIR}/test 

# SOURCE FILES FOR THE TEST EXECUTABLE (.c files needed)
EXE1_SOURCES =  ${TEST_DIR}/test.c \
				${EXTERNAL_DIR}/syng/seqhash.c \
				${EXTERNAL_DIR}/syng/utils_d.c

# OBJECT FILES IN THE BUILD DIRECTORY FROM .c FILES
EXE1_OBJECTS =  ${OBJECT_DIR}/test.o \
				${OBJECT_DIR}/seqhash.o \
				${OBJECT_DIR}/utils_d.o


PYTHON = python3
PYBIND11_INC = $(shell $(PYTHON) -m pybind11 --includes)

PYTHON_SOURCES = python/csyncmer_fast/_bindings.cpp 
PYTHON_MODULE = python/csyncmer_fast/_bindings$(shell $(PYTHON)-config --extension-suffix)

# DEFAULT TARGET FOR MAKE
all: CFLAGS += ${OPTIMIZE_FLAGS}
all : mkdir ${EXE1} # all depends from mkdir (below) and the binary ${EXE1}

#Build executable from object files
${EXE1}: ${EXE1_OBJECTS}
	${CC} ${CFLAGS} -o $@ $^ ${EXE1_LIBS}

#Compile test.c into test.o (object file)
${OBJECT_DIR}/test.o: ${TEST_DIR}/test.c
	${CC} ${CFLAGS} ${INCLUDES} -c $< -o $@

#Complie seqhash.c into seqhash.o (object file)
${OBJECT_DIR}/seqhash.o: ${EXTERNAL_DIR}/syng/seqhash.c
	${CC} ${CFLAGS} ${INCLUDES} -c $< -o $@

#Compile utils_d.c into utils_d.o (object file)
${OBJECT_DIR}/utils_d.o: ${EXTERNAL_DIR}/syng/utils_d.c
	${CC} ${CFLAGS} ${INCLUDES} -c $< -o $@

#Create necessary directories
mkdir:
	mkdir -p ${BINARY_DIR} ${OBJECT_DIR}

#Clean build artifacts
clean:
	rm -rf ${BUILD_DIR}
	rm -rf ${PYTHON_MODULE}

#Debug build
deubg: CFLAGS += ${DEBUG_FLAGS}
debug: all

#Build and run tests
test: all
	./${EXE1}

python: ${PYTHON_MODULE}
${PYTHON_MODULE}: ${PYTHON_SOURCES}
	${CC} ${CFLAGS} -I${INCLUDE_DIR} ${PYBIND11_INC} -shared -fPIC $^ -lnthash -o $@ -v
#/usr/local/lib/libnthash.a

.PHONY: all clean mkdir debug test python
