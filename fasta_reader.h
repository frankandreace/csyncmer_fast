#ifndef FASTA_READER_H
#define FASTA_READER_H

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#define INITIAL_BUF_SIZE 1024

#define STREAM_BUFSIZ (1<<12)

/* Stream type definitions */
typedef enum {
    STREAM_FILE,
    STREAM_FASTA
} stream_type;

typedef struct stream {
    stream_type type;
    int reverse;      // Not used in this simple version; reserved for future use
    long len;         // Total length of the file
    long pos;         // Current position in the file
    FILE *f;          // File pointer
    long buf_ptr;     // File offset corresponding to the start of the buffer
    size_t buf_size;  // Number of valid bytes in the buffer
    unsigned char buf[STREAM_BUFSIZ];  // Internal buffer
} stream;


/* Returns the next nucleotide from the FASTA stream:
   - Returns lowercase nucleotide letters ('a', 't', 'g', 'c', 'n').
   - Returns 0 when a new FASTA record is reached (i.e. end of a sequence).
   - Returns EOF when the file is done.
*/
static int stream_getnext_file(stream *S) {
    if (S->pos >= S->len || (S->pos - S->buf_ptr) >= (long)S->buf_size) {
        S->buf_ptr = S->pos;
        fseek(S->f, S->buf_ptr, SEEK_SET);
        S->buf_size = fread(S->buf, 1, STREAM_BUFSIZ, S->f);
        if (S->buf_size == 0)
            return EOF;
    }
    int c = S->buf[S->pos - S->buf_ptr];
    S->pos++;
    return c;
}

static int stream_getnext_fasta(stream *S) {
    int c;
    while (1) {
        c = stream_getnext_file(S);
        if (c == EOF)
            return EOF;
        if (c == '\n') {
            /* Peek at the next character */
            int next = stream_getnext_file(S);
            if (next == '>') {
                /* Rewind one character so that '>' can be processed by the next call */
                S->pos--;
                return 0;  /* End of current sequence */
            }
            continue;
        }
        if (c == '>') {
            /* Skip header line */
            while (c != '\n' && c != EOF)
                c = stream_getnext_file(S);
            continue;
        }
        /* Convert to lowercase and accept only valid nucleotide characters */
        c = tolower(c);
        if (c == 'a' || c == 't' || c == 'g' || c == 'c' || c == 'n')
            return c;
        /* Ignore other characters */
    }
}

/* Dispatch function: based on stream type call the appropriate reader. */
int stream_getnext(stream *S) {
    if (S->type == STREAM_FASTA)
        return stream_getnext_fasta(S);
    return EOF;
}

/* Opens the FASTA stream from a given FILE pointer. */
stream *stream_open_fasta(FILE *fin) {
    if (!fin)
        return NULL;
    stream *S = malloc(sizeof(stream));
    if (!S) {
        fprintf(stderr, "Not enough memory\n");
        exit(1);
    }
    S->type = STREAM_FASTA;
    S->reverse = 0;
    S->pos = 0;
    fseek(fin, 0, SEEK_END);
    S->len = ftell(fin);
    fseek(fin, 0, SEEK_SET);
    S->buf_ptr = 0;
    S->buf_size = 0;
    S->f = fin;
    return S;
}


/* reads an entire sequence from the fasta file*/
char *read_sequence(struct stream *S) {
    size_t capacity = INITIAL_BUF_SIZE;
    size_t len = 0;
    char *sequence = malloc(capacity);
    if (!sequence) {
        fprintf(stderr, "Memory allocation failure\n");
        return NULL;
    }

    int c;
    while ((c = stream_getnext(S)) != EOF && c != 0) {
        sequence[len++] = (char)c;
        if (len >= capacity) {
            capacity *= 2;
            char *tmp = realloc(sequence, capacity);
            if (!tmp) {
                free(sequence);
                fprintf(stderr, "Memory allocation failure during expansion\n");
                return NULL;
            }
            sequence = tmp;
        }
    }
    sequence[len] = '\0';
    return sequence;
}

/* Closes the FASTA stream and cleans up resources. */
void stream_close(stream *S) {
    if (S)
        free(S);
}

#endif /* FASTA_READER_H */
