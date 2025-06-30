#ifndef NTHASH_WRAPPER_H
#define NTHASH_WRAPPER_H

#include <stdint.h> // for uint64_t
#include <stddef.h>

// Checking for __uint128_t support (in GCC or CLANG)
#if defined(__SIZEOF_INT128__)
    typedef __uint128_t U128;
#else
    #error "__uint128_t is not supported by your compiler. Please use a GCC/Clang compatible compiler."
#endif

// // Forward declaration of the C++ class
// namespace nthash {
//     class NtHash;
// }
// --- C Interface Declarations ---
// These functions are declared with C linkage.
// Their *definitions* will be provided when compiled by a C++ compiler.

#ifdef __cplusplus
extern "C" { // Ensure C linkage for these functions
#endif

// // Opaque pointer to hide C++ class details from C
typedef void* NtHashHandle;

/**
 * @brief Creates a new NtHash object for a given sequence and k-mer length.
 * @param sequence The DNA/RNA sequence string.
 * @param k The k-mer length.
 * @param hash_length 1 for 64-bit canonical, 2 for 128-bit canonical.
 * @return An opaque handle to the NtHash object, or NULL on failure.
 */
NtHashHandle nthash_create(const char* sequence, size_t sequence_length, unsigned int k, unsigned int hash_length_flag);

/**
 * @brief Rolls the hash window to the next k-mer.
 * @param handle Handle to the NtHash object.
 * @return 1 if successful (more k-mers), 0 if no more k-mers.
 */
int nthash_roll(NtHashHandle handle);

int get_seq_size(NtHashHandle handle);

/**
 * @brief Gets the current forward hash (always 64-bit).
 * @param handle Handle to the NtHash object.
 * @return The 64-bit forward hash.
 */
uint64_t nthash_get_forward_hash(NtHashHandle handle);

/**
 * @brief Gets the current canonical hash.
 * Returns a 128-bit hash. If 64-bit was requested during creation,
 * the higher 64 bits of the __uint128_t will be zero.
 * @param handle Handle to the NtHash object.
 * @return The 128-bit canonical hash.
 */
U128 nthash_get_canonical_hash_128(NtHashHandle handle);

/**
 * @brief Destroys the NtHash object and frees allocated memory.
 * @param handle Opaque handle to the NtHash object.
 */
void nthash_destroy(NtHashHandle handle);

#ifdef __cplusplus
} // extern "C"
#endif

// --- C++ Implementation (Only visible to C++ compilers) ---
// This part contains the actual C++ code that implements the C interface.
// It's placed *after* the extern "C" block and is only compiled if __cplusplus is defined.

#ifdef __cplusplus

// Include necessary C++ headers here.
// These will *not* be included when a C compiler processes this header.
#include <csyncmer_fast/nthash/nthash.hpp> // The actual ntHash C++ library header
#include <vector>             // For std::vector
#include <stdexcept>          // For std::runtime_error (or other exceptions)
#include <iostream>

// Anonymous namespace or internal linkage to prevent symbol clashes
namespace{

// Helper to convert void* handle to C++ NtHash*
inline nthash::NtHash* get_nthash_obj(NtHashHandle handle){
    if (!handle) {
        throw std::runtime_error("Invalid NtHash handle.");
    }
    return static_cast<nthash::NtHash*>(handle);
}
}// anonymous namespace

// Implementations of the C functions, using the extern "C" linkage.
// Marked as 'inline' to suggest the compiler to inline them,
// which is standard for header-only libraries.
extern "C" inline NtHashHandle nthash_create(const char* sequence, size_t sequence_length, unsigned int k, unsigned int hash_length_flag){
    try{
        return new nthash::NtHash(sequence, sequence_length, hash_length_flag, k);
    } catch (const std::exception& e){
        std::cerr << "CONSTRUCTION OF NTHASH FAILED MISERABLY.\n" << std::endl;
        return nullptr;
    }
}

extern "C" inline int get_seq_size(NtHashHandle handle){
    try {
        return get_nthash_obj(handle)->get_seq_size();
    } catch (...) {return 0; }
}   

extern "C" inline int nthash_roll(NtHashHandle handle){
    try {
        return get_nthash_obj(handle)->roll() ? 1: 0;
    } catch (...) {return 0; }
}

extern "C" inline uint64_t nthash_get_forward_hash(NtHashHandle handle){
    try {
        return get_nthash_obj(handle)->get_forward_hash();
    } catch (...) {return 0; }
}

extern "C" inline U128 nthash_get_canonical_hash_128(NtHashHandle handle){
    U128 result = 0;
    try {
        const uint64_t* hashes_ptr = get_nthash_obj(handle)->hashes();
        result = ((U128)hashes_ptr[1] << 64) | hashes_ptr[0];
    } catch (...) {}
    return result;
}

extern "C" inline void nthash_destroy(NtHashHandle handle) {
    try {
        // before I deleted get_nthash_obj(handle);
        delete static_cast<nthash::NtHash*>(handle); // Use get_nthash_obj to check handle validity
    } catch (...) { /* ignore errors during destruction */ }
}

#endif // __cplusplus

#endif // NTHASH_C_API_H
