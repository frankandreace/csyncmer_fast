#pragma once

#include <cstdint>
#include <cstddef>
#include <array>
#include <cassert>
#include <vector>

namespace csyncmer_simd::hash {

/// ntHash implementation for fast rolling hash
/// Uses 32-bit hashes with R=7 bit rotation per character
namespace detail {

// Rotation constant
constexpr uint32_t R = 7;

// ntHash 32-bit lookup table (2-bit encoding: A=0, C=1, T=2, G=3)
constexpr std::array<uint32_t, 4> NTHASH_F = {
    static_cast<uint32_t>(0x3c8bfbb395c60474ULL),  // A = 0x95c60474
    static_cast<uint32_t>(0x3193c18562a02b4cULL),  // C = 0x62a02b4c
    static_cast<uint32_t>(0x20323ed082572324ULL),  // T = 0x82572324
    static_cast<uint32_t>(0x295549f54be24456ULL),  // G = 0x4be24456
};

// Complement base: A<->T (0<->2), C<->G (1<->3)
constexpr uint8_t complement(uint8_t base) {
    return base ^ 2;
}

// Complement table: c[i] = f[complement(i)]
constexpr std::array<uint32_t, 4> NTHASH_C = {
    NTHASH_F[2],  // A -> T
    NTHASH_F[3],  // C -> G
    NTHASH_F[0],  // T -> A
    NTHASH_F[1],  // G -> C
};

// Rotate left 32-bit
[[gnu::always_inline]] inline constexpr uint32_t rotl32(uint32_t x, uint32_t n) {
    n &= 31;
    return (x << n) | (x >> (32 - n));
}

// Rotate right 32-bit
[[gnu::always_inline]] inline constexpr uint32_t rotr32(uint32_t x, uint32_t n) {
    n &= 31;
    return (x >> n) | (x << (32 - n));
}

// Fixed rotation by R=7 bits (hot path)
[[gnu::always_inline]] inline uint32_t rotl7(uint32_t x) {
    return (x << 7) | (x >> 25);
}

// Pre-rotated tables (rotated by (k-1)*R)
inline std::array<uint32_t, 4> make_f_rot(size_t k) {
    uint32_t rot = static_cast<uint32_t>((k - 1) * R);
    return {
        rotl32(NTHASH_F[0], rot),
        rotl32(NTHASH_F[1], rot),
        rotl32(NTHASH_F[2], rot),
        rotl32(NTHASH_F[3], rot),
    };
}

inline std::array<uint32_t, 4> make_c_rot(size_t k) {
    uint32_t rot = static_cast<uint32_t>((k - 1) * R);
    return {
        rotl32(NTHASH_C[0], rot),
        rotl32(NTHASH_C[1], rot),
        rotl32(NTHASH_C[2], rot),
        rotl32(NTHASH_C[3], rot),
    };
}

// Initial value of hashing k-1 zeros
inline uint32_t make_fw_init(size_t k) {
    uint32_t fw = 0;
    for (size_t i = 0; i < k - 1; ++i) {
        fw = rotl32(fw, R) ^ NTHASH_F[0];
    }
    return fw;
}

inline uint32_t make_rc_init(size_t k, const std::array<uint32_t, 4>& c_rot) {
    uint32_t rc = 0;
    for (size_t i = 0; i < k - 1; ++i) {
        rc = rotr32(rc, R) ^ c_rot[0];
    }
    return rc;
}

} // namespace detail

/// Forward-only ntHash for DNA sequences
/// 32-bit hash with R=7 rotation
class NtHash {
public:
    NtHash() = default;

    explicit NtHash(size_t k)
        : k_(k)
        , f_rot_(detail::make_f_rot(k))
        , fw_init_(detail::make_fw_init(k))
        , fw_(fw_init_)
    {}

    /// Initialize and return first hash
    /// CORRECT APPROACH: Compute hash directly, then set up state for rolling
    template<typename SeqT>
    uint32_t init(const SeqT& seq, size_t start) {
        // Compute hash directly: H = rotl^(k-1)(F[s[0]]) ^ rotl^(k-2)(F[s[1]]) ^ ... ^ F[s[k-1]]
        uint32_t hash = 0;
        for (size_t i = 0; i < k_; ++i) {
            hash = detail::rotl32(hash, detail::R) ^ detail::NTHASH_F[seq[start + i]];
        }

        // Set up state for rolling: fw_ = hash ^ f_rot[first_base]
        fw_ = hash ^ f_rot_[seq[start]];

        return hash;
    }

    /// Roll: remove out_base, add in_base, new_first_base is the first base of new s-mer
    /// After rolling from position i to i+1:
    ///   - out_base = seq[i] (being removed)
    ///   - in_base = seq[i + k] (being added)
    ///   - new_first_base = seq[i + 1] (first base of new s-mer, needed for state)
    [[gnu::always_inline]] uint32_t roll(uint8_t out_base, uint8_t in_base, uint8_t new_first_base) {
        (void)out_base;  // Not needed - it's already encoded in fw_ state
        uint32_t fw_out = detail::rotl7(fw_) ^ detail::NTHASH_F[in_base];
        fw_ = fw_out ^ f_rot_[new_first_base];
        return fw_out;
    }

    /// Legacy roll API (2 args) - only works when new_first_base == out_base
    /// This is incorrect for general use! Kept for compatibility.
    [[gnu::always_inline, deprecated("Use roll(out, in, new_first) instead")]]
    uint32_t roll_legacy(uint8_t out_base, uint8_t in_base) {
        uint32_t fw_out = detail::rotl7(fw_) ^ detail::NTHASH_F[in_base];
        fw_ = fw_out ^ f_rot_[out_base];
        return fw_out;
    }

    [[nodiscard]] uint32_t hash() const { return fw_; }

private:
    size_t k_ = 0;
    std::array<uint32_t, 4> f_rot_{};
    uint32_t fw_init_ = 0;
    uint32_t fw_ = 0;
};

/// Canonical ntHash - same hash for k-mer and its reverse complement
class CanonicalNtHash {
public:
    CanonicalNtHash() = default;

    explicit CanonicalNtHash(size_t k)
        : k_(k)
        , f_rot_(detail::make_f_rot(k))
        , c_rot_(detail::make_c_rot(k))
        , fw_init_(detail::make_fw_init(k))
        , rc_init_(detail::make_rc_init(k, c_rot_))
        , fw_(fw_init_)
        , rc_(rc_init_)
    {}

    /// Initialize and return first canonical hash
    /// CORRECT APPROACH: Compute hash directly, then set up state for rolling
    template<typename SeqT>
    uint32_t init(const SeqT& seq, size_t start) {
        // Compute forward hash directly
        uint32_t fw_hash = 0;
        for (size_t i = 0; i < k_; ++i) {
            fw_hash = detail::rotl32(fw_hash, detail::R) ^ detail::NTHASH_F[seq[start + i]];
        }

        // Compute reverse complement hash directly
        uint32_t rc_hash = 0;
        for (size_t i = 0; i < k_; ++i) {
            rc_hash = detail::rotr32(rc_hash, detail::R) ^ c_rot_[seq[start + i]];
        }

        // Set up state for rolling
        uint8_t first_base = seq[start];
        fw_ = fw_hash ^ f_rot_[first_base];
        rc_ = rc_hash ^ detail::NTHASH_C[first_base];

        return fw_hash + rc_hash;
    }

    uint32_t roll(uint8_t out_base, uint8_t in_base, uint8_t new_first_base) {
        (void)out_base;  // Not needed - already encoded in state
        uint32_t fw_out = detail::rotl32(fw_, detail::R) ^ detail::NTHASH_F[in_base];
        fw_ = fw_out ^ f_rot_[new_first_base];
        uint32_t rc_out = detail::rotr32(rc_, detail::R) ^ c_rot_[in_base];
        rc_ = rc_out ^ detail::NTHASH_C[new_first_base];
        return fw_out + rc_out;
    }

    [[nodiscard]] uint32_t hash() const { return fw_ + rc_; }

private:
    size_t k_ = 0;
    std::array<uint32_t, 4> f_rot_{};
    std::array<uint32_t, 4> c_rot_{};
    uint32_t fw_init_ = 0;
    uint32_t rc_init_ = 0;
    uint32_t fw_ = 0;
    uint32_t rc_ = 0;
};

} // namespace csyncmer_simd::hash
