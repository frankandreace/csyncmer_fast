#pragma once

#include <array>
#include <cstdint>
#include <algorithm>
#include <limits>

// Platform detection
#if defined(__ARM_NEON) || defined(__ARM_NEON__)
    #include <arm_neon.h>
    #define CSYNCMER_SIMD_NEON 1
#elif defined(__AVX2__)
    #include <immintrin.h>
    #define CSYNCMER_SIMD_AVX2 1
#endif

namespace csyncmer_simd {

constexpr size_t LANES = 8;

/// SIMD vector of 8 x uint32_t
/// On ARM NEON, uses 2x 128-bit registers
/// On AVX2, uses 1x 256-bit register
/// Otherwise, scalar fallback with std::array
#if defined(CSYNCMER_SIMD_NEON)

struct u32x8 {
    uint32x4_t lo;
    uint32x4_t hi;

    u32x8() = default;
    u32x8(uint32x4_t l, uint32x4_t h) : lo(l), hi(h) {}

    static u32x8 splat(uint32_t v) {
        return {vdupq_n_u32(v), vdupq_n_u32(v)};
    }

    static u32x8 zero() {
        return {vdupq_n_u32(0), vdupq_n_u32(0)};
    }

    static u32x8 max_value() {
        return splat(std::numeric_limits<uint32_t>::max());
    }

    // Load from array
    static u32x8 load(const uint32_t* ptr) {
        return {vld1q_u32(ptr), vld1q_u32(ptr + 4)};
    }

    static u32x8 from_array(const std::array<uint32_t, 8>& arr) {
        return load(arr.data());
    }

    template<typename F>
    static u32x8 from_fn(F&& f) {
        std::array<uint32_t, 8> arr;
        for (size_t i = 0; i < 8; ++i) arr[i] = f(i);
        return from_array(arr);
    }

    // Store to array
    void store(uint32_t* ptr) const {
        vst1q_u32(ptr, lo);
        vst1q_u32(ptr + 4, hi);
    }

    // Element access (for debugging/testing)
    uint32_t operator[](size_t i) const {
        uint32_t arr[8];
        store(arr);
        return arr[i];
    }

    uint32_t get(size_t i) const { return (*this)[i]; }

    // Min
    u32x8 min(const u32x8& other) const {
        return {vminq_u32(lo, other.lo), vminq_u32(hi, other.hi)};
    }

    // Max
    u32x8 max(const u32x8& other) const {
        return {vmaxq_u32(lo, other.lo), vmaxq_u32(hi, other.hi)};
    }

    // Add
    u32x8 operator+(const u32x8& other) const {
        return {vaddq_u32(lo, other.lo), vaddq_u32(hi, other.hi)};
    }

    // Subtract
    u32x8 operator-(const u32x8& other) const {
        return {vsubq_u32(lo, other.lo), vsubq_u32(hi, other.hi)};
    }

    u32x8& operator+=(const u32x8& other) {
        lo = vaddq_u32(lo, other.lo);
        hi = vaddq_u32(hi, other.hi);
        return *this;
    }

    u32x8& operator-=(const u32x8& other) {
        lo = vsubq_u32(lo, other.lo);
        hi = vsubq_u32(hi, other.hi);
        return *this;
    }

    // Bitwise AND
    u32x8 operator&(const u32x8& other) const {
        return {vandq_u32(lo, other.lo), vandq_u32(hi, other.hi)};
    }

    // Bitwise OR
    u32x8 operator|(const u32x8& other) const {
        return {vorrq_u32(lo, other.lo), vorrq_u32(hi, other.hi)};
    }

    // Bitwise XOR
    u32x8 operator^(const u32x8& other) const {
        return {veorq_u32(lo, other.lo), veorq_u32(hi, other.hi)};
    }

    // Bitwise NOT
    u32x8 operator~() const {
        return {vmvnq_u32(lo), vmvnq_u32(hi)};
    }

    // Shift left by scalar
    u32x8 shl(int shift) const {
        int32x4_t shift_vec = vdupq_n_s32(shift);
        return {
            vshlq_u32(lo, shift_vec),
            vshlq_u32(hi, shift_vec)
        };
    }

    // Compare equal (returns mask)
    u32x8 cmp_eq(const u32x8& other) const {
        return {vceqq_u32(lo, other.lo), vceqq_u32(hi, other.hi)};
    }

    // Compare greater than (returns mask)
    u32x8 cmp_gt(const u32x8& other) const {
        return {vcgtq_u32(lo, other.lo), vcgtq_u32(hi, other.hi)};
    }

    // Blend: select from a where mask is true, else b
    static u32x8 blend(const u32x8& mask, const u32x8& a, const u32x8& b) {
        return {
            vbslq_u32(mask.lo, a.lo, b.lo),
            vbslq_u32(mask.hi, a.hi, b.hi)
        };
    }

    // Extract 8-bit mask from comparison results (1 bit per lane)
    uint8_t movemask() const {
        // NEON doesn't have movemask, emulate it
        uint32x4_t shift_lo = vshrq_n_u32(lo, 31);
        uint32x4_t shift_hi = vshrq_n_u32(hi, 31);

        uint64_t lo_bits = vgetq_lane_u32(shift_lo, 0)
                        | (vgetq_lane_u32(shift_lo, 1) << 1)
                        | (vgetq_lane_u32(shift_lo, 2) << 2)
                        | (vgetq_lane_u32(shift_lo, 3) << 3);
        uint64_t hi_bits = vgetq_lane_u32(shift_hi, 0)
                        | (vgetq_lane_u32(shift_hi, 1) << 1)
                        | (vgetq_lane_u32(shift_hi, 2) << 2)
                        | (vgetq_lane_u32(shift_hi, 3) << 3);

        return static_cast<uint8_t>(lo_bits | (hi_bits << 4));
    }

    // Horizontal reduce (get minimum across all lanes)
    uint32_t reduce_min() const {
        uint32x4_t m = vminq_u32(lo, hi);
        m = vpminq_u32(m, m);
        m = vpminq_u32(m, m);
        return vgetq_lane_u32(m, 0);
    }

    // To array
    std::array<uint32_t, 8> to_array() const {
        std::array<uint32_t, 8> arr;
        store(arr.data());
        return arr;
    }
};

inline u32x8 U32X8_ZERO() { return u32x8::zero(); }
inline u32x8 U32X8_MAX() { return u32x8::max_value(); }
inline u32x8 U32X8_ONE() { return u32x8::splat(1); }

#elif defined(CSYNCMER_SIMD_AVX2)

struct u32x8 {
    __m256i v;

    u32x8() = default;
    explicit u32x8(__m256i val) : v(val) {}

    static u32x8 splat(uint32_t val) {
        return u32x8{_mm256_set1_epi32(static_cast<int32_t>(val))};
    }

    static u32x8 zero() {
        return u32x8{_mm256_setzero_si256()};
    }

    static u32x8 max_value() {
        return splat(std::numeric_limits<uint32_t>::max());
    }

    static u32x8 load(const uint32_t* ptr) {
        return u32x8{_mm256_loadu_si256(reinterpret_cast<const __m256i*>(ptr))};
    }

    static u32x8 from_array(const std::array<uint32_t, 8>& arr) {
        return load(arr.data());
    }

    template<typename F>
    static u32x8 from_fn(F&& f) {
        std::array<uint32_t, 8> arr;
        for (size_t i = 0; i < 8; ++i) arr[i] = f(i);
        return from_array(arr);
    }

    void store(uint32_t* ptr) const {
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(ptr), v);
    }

    uint32_t operator[](size_t i) const {
        alignas(32) uint32_t arr[8];
        store(arr);
        return arr[i];
    }

    uint32_t get(size_t i) const { return (*this)[i]; }

    u32x8 min(const u32x8& other) const {
        return u32x8{_mm256_min_epu32(v, other.v)};
    }

    u32x8 max(const u32x8& other) const {
        return u32x8{_mm256_max_epu32(v, other.v)};
    }

    u32x8 operator+(const u32x8& other) const {
        return u32x8{_mm256_add_epi32(v, other.v)};
    }

    u32x8 operator-(const u32x8& other) const {
        return u32x8{_mm256_sub_epi32(v, other.v)};
    }

    u32x8& operator+=(const u32x8& other) {
        v = _mm256_add_epi32(v, other.v);
        return *this;
    }

    u32x8& operator-=(const u32x8& other) {
        v = _mm256_sub_epi32(v, other.v);
        return *this;
    }

    u32x8 operator&(const u32x8& other) const {
        return u32x8{_mm256_and_si256(v, other.v)};
    }

    u32x8 operator|(const u32x8& other) const {
        return u32x8{_mm256_or_si256(v, other.v)};
    }

    u32x8 operator^(const u32x8& other) const {
        return u32x8{_mm256_xor_si256(v, other.v)};
    }

    u32x8 operator~() const {
        return u32x8{_mm256_xor_si256(v, _mm256_set1_epi32(-1))};
    }

    u32x8 shl(int shift) const {
        return u32x8{_mm256_slli_epi32(v, shift)};
    }

    u32x8 cmp_eq(const u32x8& other) const {
        return u32x8{_mm256_cmpeq_epi32(v, other.v)};
    }

    u32x8 cmp_gt(const u32x8& other) const {
        // AVX2 only has signed comparison, need to flip sign bit
        __m256i sign = _mm256_set1_epi32(static_cast<int32_t>(0x80000000u));
        __m256i a = _mm256_xor_si256(v, sign);
        __m256i b = _mm256_xor_si256(other.v, sign);
        return u32x8{_mm256_cmpgt_epi32(a, b)};
    }

    static u32x8 blend(const u32x8& mask, const u32x8& a, const u32x8& b) {
        return u32x8{_mm256_blendv_epi8(b.v, a.v, mask.v)};
    }

    uint8_t movemask() const {
        return static_cast<uint8_t>(_mm256_movemask_ps(_mm256_castsi256_ps(v)));
    }

    uint32_t reduce_min() const {
        // Get upper and lower 128-bit halves
        __m128i lo = _mm256_castsi256_si128(v);
        __m128i hi = _mm256_extracti128_si256(v, 1);
        __m128i m = _mm_min_epu32(lo, hi);
        m = _mm_min_epu32(m, _mm_shuffle_epi32(m, 0x4E)); // 2,3,0,1
        m = _mm_min_epu32(m, _mm_shuffle_epi32(m, 0xB1)); // 1,0,3,2
        return static_cast<uint32_t>(_mm_cvtsi128_si32(m));
    }

    std::array<uint32_t, 8> to_array() const {
        alignas(32) std::array<uint32_t, 8> arr;
        store(arr.data());
        return arr;
    }
};

inline u32x8 U32X8_ZERO() { return u32x8::zero(); }
inline u32x8 U32X8_MAX() { return u32x8::max_value(); }
inline u32x8 U32X8_ONE() { return u32x8::splat(1); }

#else
// Scalar fallback

struct u32x8 {
    std::array<uint32_t, 8> data;

    u32x8() = default;
    explicit u32x8(const std::array<uint32_t, 8>& d) : data(d) {}

    static u32x8 splat(uint32_t v) {
        return u32x8{{v, v, v, v, v, v, v, v}};
    }

    static u32x8 zero() {
        return splat(0);
    }

    static u32x8 max_value() {
        return splat(std::numeric_limits<uint32_t>::max());
    }

    static u32x8 load(const uint32_t* ptr) {
        u32x8 result;
        for (size_t i = 0; i < 8; ++i) result.data[i] = ptr[i];
        return result;
    }

    static u32x8 from_array(const std::array<uint32_t, 8>& arr) {
        return u32x8{arr};
    }

    template<typename F>
    static u32x8 from_fn(F&& f) {
        std::array<uint32_t, 8> arr;
        for (size_t i = 0; i < 8; ++i) arr[i] = f(i);
        return from_array(arr);
    }

    void store(uint32_t* ptr) const {
        for (size_t i = 0; i < 8; ++i) ptr[i] = data[i];
    }

    uint32_t operator[](size_t i) const { return data[i]; }
    uint32_t& operator[](size_t i) { return data[i]; }
    uint32_t get(size_t i) const { return data[i]; }

    u32x8 min(const u32x8& other) const {
        u32x8 r;
        for (size_t i = 0; i < 8; ++i) r.data[i] = std::min(data[i], other.data[i]);
        return r;
    }

    u32x8 max(const u32x8& other) const {
        u32x8 r;
        for (size_t i = 0; i < 8; ++i) r.data[i] = std::max(data[i], other.data[i]);
        return r;
    }

    u32x8 operator+(const u32x8& other) const {
        u32x8 r;
        for (size_t i = 0; i < 8; ++i) r.data[i] = data[i] + other.data[i];
        return r;
    }

    u32x8 operator-(const u32x8& other) const {
        u32x8 r;
        for (size_t i = 0; i < 8; ++i) r.data[i] = data[i] - other.data[i];
        return r;
    }

    u32x8& operator+=(const u32x8& other) {
        for (size_t i = 0; i < 8; ++i) data[i] += other.data[i];
        return *this;
    }

    u32x8& operator-=(const u32x8& other) {
        for (size_t i = 0; i < 8; ++i) data[i] -= other.data[i];
        return *this;
    }

    u32x8 operator&(const u32x8& other) const {
        u32x8 r;
        for (size_t i = 0; i < 8; ++i) r.data[i] = data[i] & other.data[i];
        return r;
    }

    u32x8 operator|(const u32x8& other) const {
        u32x8 r;
        for (size_t i = 0; i < 8; ++i) r.data[i] = data[i] | other.data[i];
        return r;
    }

    u32x8 operator^(const u32x8& other) const {
        u32x8 r;
        for (size_t i = 0; i < 8; ++i) r.data[i] = data[i] ^ other.data[i];
        return r;
    }

    u32x8 operator~() const {
        u32x8 r;
        for (size_t i = 0; i < 8; ++i) r.data[i] = ~data[i];
        return r;
    }

    u32x8 shl(int shift) const {
        u32x8 r;
        for (size_t i = 0; i < 8; ++i) r.data[i] = data[i] << shift;
        return r;
    }

    u32x8 cmp_eq(const u32x8& other) const {
        u32x8 r;
        for (size_t i = 0; i < 8; ++i) {
            r.data[i] = data[i] == other.data[i] ? ~0u : 0u;
        }
        return r;
    }

    u32x8 cmp_gt(const u32x8& other) const {
        u32x8 r;
        for (size_t i = 0; i < 8; ++i) {
            r.data[i] = data[i] > other.data[i] ? ~0u : 0u;
        }
        return r;
    }

    static u32x8 blend(const u32x8& mask, const u32x8& a, const u32x8& b) {
        u32x8 r;
        for (size_t i = 0; i < 8; ++i) {
            r.data[i] = mask.data[i] ? a.data[i] : b.data[i];
        }
        return r;
    }

    uint8_t movemask() const {
        uint8_t m = 0;
        for (size_t i = 0; i < 8; ++i) {
            if (data[i] & 0x80000000) m |= (1u << i);
        }
        return m;
    }

    uint32_t reduce_min() const {
        return *std::min_element(data.begin(), data.end());
    }

    std::array<uint32_t, 8> to_array() const {
        return data;
    }
};

inline u32x8 U32X8_ZERO() { return u32x8::zero(); }
inline u32x8 U32X8_MAX() { return u32x8::max_value(); }
inline u32x8 U32X8_ONE() { return u32x8::splat(1); }

#endif

/// Signed version for comparison operations
struct i32x8 {
#if defined(CSYNCMER_SIMD_NEON)
    int32x4_t lo;
    int32x4_t hi;

    i32x8() = default;
    i32x8(int32x4_t l, int32x4_t h) : lo(l), hi(h) {}

    static i32x8 splat(int32_t v) {
        return {vdupq_n_s32(v), vdupq_n_s32(v)};
    }

    static i32x8 zero() { return splat(0); }

    i32x8 operator+(const i32x8& other) const {
        return {vaddq_s32(lo, other.lo), vaddq_s32(hi, other.hi)};
    }

    i32x8 operator-(const i32x8& other) const {
        return {vsubq_s32(lo, other.lo), vsubq_s32(hi, other.hi)};
    }

    i32x8 operator&(const i32x8& other) const {
        return {vandq_s32(lo, other.lo), vandq_s32(hi, other.hi)};
    }

    u32x8 cmp_gt(const i32x8& other) const {
        return {
            vreinterpretq_u32_s32(vcgtq_s32(lo, other.lo)),
            vreinterpretq_u32_s32(vcgtq_s32(hi, other.hi))
        };
    }

    int32_t operator[](size_t i) const {
        int32_t arr[8];
        vst1q_s32(arr, lo);
        vst1q_s32(arr + 4, hi);
        return arr[i];
    }
#elif defined(CSYNCMER_SIMD_AVX2)
    __m256i v;

    i32x8() = default;
    explicit i32x8(__m256i val) : v(val) {}

    static i32x8 splat(int32_t val) {
        return i32x8{_mm256_set1_epi32(val)};
    }

    static i32x8 zero() { return i32x8{_mm256_setzero_si256()}; }

    i32x8 operator+(const i32x8& other) const {
        return i32x8{_mm256_add_epi32(v, other.v)};
    }

    i32x8 operator-(const i32x8& other) const {
        return i32x8{_mm256_sub_epi32(v, other.v)};
    }

    i32x8 operator&(const i32x8& other) const {
        return i32x8{_mm256_and_si256(v, other.v)};
    }

    u32x8 cmp_gt(const i32x8& other) const {
        return u32x8{_mm256_cmpgt_epi32(v, other.v)};
    }

    int32_t operator[](size_t i) const {
        alignas(32) int32_t arr[8];
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(arr), v);
        return arr[i];
    }
#else
    std::array<int32_t, 8> data;

    i32x8() = default;
    explicit i32x8(const std::array<int32_t, 8>& d) : data(d) {}

    static i32x8 splat(int32_t v) {
        return i32x8{{v, v, v, v, v, v, v, v}};
    }

    static i32x8 zero() { return splat(0); }

    i32x8 operator+(const i32x8& other) const {
        i32x8 r;
        for (size_t i = 0; i < 8; ++i) r.data[i] = data[i] + other.data[i];
        return r;
    }

    i32x8 operator-(const i32x8& other) const {
        i32x8 r;
        for (size_t i = 0; i < 8; ++i) r.data[i] = data[i] - other.data[i];
        return r;
    }

    i32x8 operator&(const i32x8& other) const {
        i32x8 r;
        for (size_t i = 0; i < 8; ++i) r.data[i] = data[i] & other.data[i];
        return r;
    }

    u32x8 cmp_gt(const i32x8& other) const {
        u32x8 r;
        for (size_t i = 0; i < 8; ++i) {
            r.data[i] = data[i] > other.data[i] ? ~0u : 0u;
        }
        return r;
    }

    int32_t operator[](size_t i) const { return data[i]; }
    int32_t& operator[](size_t i) { return data[i]; }
#endif
};

// Type cast between u32x8 and i32x8
inline i32x8 reinterpret_i32(const u32x8& v) {
#if defined(CSYNCMER_SIMD_NEON)
    return {vreinterpretq_s32_u32(v.lo), vreinterpretq_s32_u32(v.hi)};
#elif defined(CSYNCMER_SIMD_AVX2)
    return i32x8{v.v};
#else
    i32x8 r;
    for (size_t i = 0; i < 8; ++i) {
        r.data[i] = static_cast<int32_t>(v.data[i]);
    }
    return r;
#endif
}

inline u32x8 reinterpret_u32(const i32x8& v) {
#if defined(CSYNCMER_SIMD_NEON)
    return {vreinterpretq_u32_s32(v.lo), vreinterpretq_u32_s32(v.hi)};
#elif defined(CSYNCMER_SIMD_AVX2)
    return u32x8{v.v};
#else
    u32x8 r;
    for (size_t i = 0; i < 8; ++i) {
        r.data[i] = static_cast<uint32_t>(v.data[i]);
    }
    return r;
#endif
}

/// Validity masks for masking out invalid elements in a u32x8 vector
/// VALIDITY_MASKS[n] has first n elements as 0, rest as 0xFFFFFFFF
/// Used to invalidate elements beyond a certain count using blend
namespace detail {

inline std::array<u32x8, 9> generate_validity_masks() {
    std::array<u32x8, 9> masks;
    for (size_t n = 0; n <= 8; ++n) {
        std::array<uint32_t, 8> arr;
        for (size_t i = 0; i < 8; ++i) {
            arr[i] = (i < n) ? 0 : UINT32_MAX;
        }
        masks[n] = u32x8::from_array(arr);
    }
    return masks;
}

} // namespace detail

/// Pre-computed validity masks
/// VALIDITY_MASKS[n]: first n elements = 0 (keep), rest = MAX (invalidate)
inline const std::array<u32x8, 9>& get_validity_masks() {
    static const auto masks = detail::generate_validity_masks();
    return masks;
}

} // namespace csyncmer_simd
