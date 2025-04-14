CC = gcc
CFLAGS = -march=native -O3
DEBUG_FLAGS = -DDEBUG -g
OUTDIR = ./bin


EXE1 = bin/test

EXE1_SRCS = test.c 

EXE1_OBJS = $(EXE1_SRCS:.c=.o)

all: mkdir
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