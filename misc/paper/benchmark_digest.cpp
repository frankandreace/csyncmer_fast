// Standalone digest syncmer benchmark for FASTA/FASTQ files.
// Compile: g++ -std=c++17 -O3 -o benchmark_digest benchmark_digest.cpp
//          -I <digest>/build/include -L <digest>/build/lib -lnthash
// Usage:   ./benchmark_digest <input> <k_small> <large_window> [runs]
//
// Prints one line: "throughput_gbps <value>"

#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>

#include <digest/syncmer.hpp>
#include <digest/data_structure.hpp>

// Simple FASTA/FASTQ parser: returns all sequences concatenated into `seqs`.
static void read_sequences(const char *path, std::vector<std::string> &seqs) {
    std::ifstream in(path);
    if (!in) { fprintf(stderr, "Cannot open %s\n", path); exit(1); }

    std::string line, seq;
    bool is_fastq = false;
    enum { HEADER, SEQ, PLUS, QUAL } state = HEADER;

    while (std::getline(in, line)) {
        if (line.empty()) continue;

        if (state == HEADER) {
            if (line[0] == '>') {
                if (!seq.empty()) { seqs.push_back(std::move(seq)); seq.clear(); }
                state = SEQ;
            } else if (line[0] == '@') {
                if (!seq.empty()) { seqs.push_back(std::move(seq)); seq.clear(); }
                is_fastq = true;
                state = SEQ;
            }
        } else if (state == SEQ) {
            if (line[0] == '>') {
                if (!seq.empty()) { seqs.push_back(std::move(seq)); seq.clear(); }
                state = SEQ; // next sequence header (FASTA)
            } else if (line[0] == '+' && is_fastq) {
                if (!seq.empty()) { seqs.push_back(std::move(seq)); seq.clear(); }
                state = QUAL;
            } else {
                seq += line;
            }
        } else if (state == QUAL) {
            // skip quality line
            state = HEADER;
        }
    }
    if (!seq.empty()) seqs.push_back(std::move(seq));
}

int main(int argc, char **argv) {
    if (argc < 4) {
        fprintf(stderr, "Usage: %s <input> <k_small> <large_window> [runs]\n", argv[0]);
        return 1;
    }

    const char *input = argv[1];
    unsigned k_small = atoi(argv[2]);
    unsigned large_window = atoi(argv[3]);
    int runs = argc > 4 ? atoi(argv[4]) : 3;

    fprintf(stderr, "Reading %s ...\n", input);
    std::vector<std::string> seqs;
    read_sequences(input, seqs);

    size_t total_bases = 0;
    for (auto &s : seqs) total_bases += s.size();
    fprintf(stderr, "Read %zu sequences, %zu bases (%.2f GB)\n",
            seqs.size(), total_bases, total_bases / 1e9);

    // Warm-up + benchmark
    double best_secs = 1e30;
    size_t total_syncmers = 0;

    for (int r = 0; r < runs; r++) {
        std::vector<uint32_t> vec;
        vec.reserve(total_bases / 10); // rough estimate

        auto t0 = std::chrono::high_resolution_clock::now();

        for (auto &s : seqs) {
            if (s.size() < k_small + large_window) continue;
            digest::Syncmer<digest::BadCharPolicy::SKIPOVER, digest::ds::Adaptive>
                dig(s, k_small, large_window, 0, digest::MinimizedHashType::CANON);
            dig.roll_minimizer(s.size(), vec);
        }

        auto t1 = std::chrono::high_resolution_clock::now();
        double secs = std::chrono::duration<double>(t1 - t0).count();

        fprintf(stderr, "  Run %d: %.3f s  (%zu syncmers)\n", r + 1, secs, vec.size());

        if (secs < best_secs) {
            best_secs = secs;
            total_syncmers = vec.size();
        }
    }

    double gbps = (total_bases / 1e9) / best_secs;
    fprintf(stderr, "Best: %.3f s => %.4f GB/s  (%zu syncmers)\n",
            best_secs, gbps, total_syncmers);

    printf("%.4f\n", gbps);
    return 0;
}
