// test_simple.cpp - enhanced version
#include <iostream>
#include <cstring>
#include "csyncmer_fast/nthash/nthash.hpp"

int main() {
    const char* seq = "AACCAAGGTTAA";
    size_t len = strlen(seq);
    
    // Try different combinations
    struct TestCase {
        unsigned num_hashes;
        unsigned k;
        const char* desc;
    };
    
    TestCase cases[] = {
        {1, 9, "1 hash, k=9"},
        {1, 8, "1 hash, k=8"},
        {1, 7, "1 hash, k=7"},
        {1, 6, "1 hash, k=6"},
        {2, 5, "2 hashes, k=5"},
        {3, 4, "3 hashes, k=4"},
        {3, 3, "3 hashes, k=3"},
        {3, 2, "3 hashes, k=2"},
    };
    
    for (const auto& tc : cases) {
        std::cout << "\nTesting: " << tc.desc << std::endl;
        try {
            nthash::NtHash hasher(seq, len, tc.num_hashes, tc.k);
            std::cout << "  Created OK, seq_size=" << hasher.get_seq_size() << std::endl;
            
            bool result = hasher.roll();
            std::cout << "  Roll result: " << result << std::endl;
            
            if (result) {
                std::cout << "  Success! pos=" << hasher.get_pos() << std::endl;
            }
        } catch (const std::exception& e) {
            std::cerr << "  Exception: " << e.what() << std::endl;
        } catch (...) {
            std::cerr << "  Unknown exception" << std::endl;
        }
    }
    
    return 0;
}