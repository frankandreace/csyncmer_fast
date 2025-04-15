CC = gcc
CFLAGS = -march=native -O3 -fno-tree-vectorize
DEBUG_FLAGS = -DDEBUG -g
OUTDIR = ./bin

EXE1_LIBS = -lz

EXE1 = bin/test

EXE1_SRCS = test.c syng/seqhash.c syng/utils_d.c

EXE1_OBJS = $(EXE1_SRCS:.c=.o)

all: mkdir
all: clean
all: $(EXE1)

$(EXE1): $(EXE1_OBJS)
	$(CC) $(CFLAGS) -o $(EXE1) $(EXE1_OBJS) $(EXE1_LIBS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

mkdir:
	mkdir -p $(OUTDIR)

clean:
	rm -f $(EXE1) $(EXE1_OBJS)

.PHONY: all clean