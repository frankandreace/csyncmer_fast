#pragma once

#include "vec.hpp"
#include <vector>
#include <cassert>
#include <cstddef>

namespace csyncmer_simd {

/// SIMD Ring buffer for sliding window
class RingBufSimd {
public:
    RingBufSimd() = default;

    void init(size_t w, u32x8 fill_value) {
        w_ = w;
        idx_ = 0;
        data_.resize(w, fill_value);
    }

    void push(u32x8 val) {
        data_[idx_] = val;
        idx_ = (idx_ + 1) % w_;
    }

    size_t idx() const { return idx_; }
    size_t size() const { return w_; }

    u32x8& operator[](size_t i) { return data_[i]; }
    const u32x8& operator[](size_t i) const { return data_[i]; }

    u32x8* begin() { return data_.data(); }
    u32x8* end() { return data_.data() + w_; }

private:
    size_t w_ = 0;
    size_t idx_ = 0;
    std::vector<u32x8> data_;
};

/// SIMD Sliding Window Minimum using Two-Stacks algorithm
/// Processes 8 parallel windows (one per SIMD lane)
///
/// Template parameter LEFT:
///   true  = leftmost minimum (ties broken by smallest position)
///   false = rightmost minimum (ties broken by largest position)
template<bool LEFT = true>
class SlidingMinSimd {
public:
    SlidingMinSimd() = default;

    /// Initialize for window size w over total_len elements per chunk
    void init(size_t w, size_t total_len) {
        assert(w > 0);
        assert(w < (1u << 15));  // Position fits in 16 bits

        w_ = w;
        total_len_ = total_len;

        // Pre-computed masks
        val_mask_ = u32x8::splat(0xFFFF0000);
        pos_mask_ = u32x8::splat(0x0000FFFF);
        max_pos_ = u32x8::splat((1u << 16) - 1);

        // Position counter (starts at 0)
        pos_ = U32X8_ZERO();

        // Position offset per lane: each lane tracks its own internal position
        // We don't need external offsets since get_window_offset works within each lane
        pos_offset_ = U32X8_ZERO();

        // Delta for position reset
        delta_ = u32x8::splat((1u << 16) - static_cast<uint32_t>(w));

        // Ring buffer for suffix computation
        ring_buf_.init(w, U32X8_MAX());

        // Prefix minimum starts at infinity
        prefix_min_ = U32X8_MAX();

        // Initialize ring buffer with k-1 "zeros" (position 0, value infinity)
        fw_init_ = U32X8_ZERO();
        for (size_t i = 0; i < w - 1; ++i) {
            u32x8 elem = val_mask_;  // High value with position 0
            if constexpr (!LEFT) {
                elem = ~elem & val_mask_;  // Invert for RIGHT
            }
            fw_init_ = simd_min(fw_init_, elem);
        }
    }

    /// Push a new hash value (upper 16 bits) and return window minimum position
    u32x8 push(u32x8 hash_value) {
        // Check for position overflow
        auto overflow = pos_.cmp_eq(max_pos_);
        if (overflow.movemask() != 0) {
            reset_positions();
        }

        // Encode: upper 16 bits = hash, lower 16 bits = position
        u32x8 elem = (hash_value & val_mask_) | pos_;

        // For RIGHT tie-breaking, invert the value bits
        if constexpr (!LEFT) {
            elem = (~hash_value & val_mask_) | pos_;
        }

        // Increment position
        pos_ += U32X8_ONE();

        // Push to ring buffer
        ring_buf_.push(elem);

        // Update prefix minimum
        prefix_min_ = simd_min(prefix_min_, elem);

        // When ring buffer wraps, compute suffix minima
        if (ring_buf_.idx() == 0) {
            compute_suffix_minima(elem);
        }

        // Get suffix minimum from current position
        u32x8 suffix_min = ring_buf_[ring_buf_.idx()];

        // Window minimum = min(prefix, suffix)
        u32x8 window_min = simd_min(prefix_min_, suffix_min);

        // Extract position and add offset
        return (window_min & pos_mask_) + pos_offset_;
    }

    /// Get the position within the window (0 = oldest, w-1 = newest)
    /// This is needed for closed syncmer detection
    u32x8 get_window_offset(u32x8 min_pos) const {
        // Current position is the position of the newest element
        // min_pos includes pos_offset_ (returned by push as: raw_pos + pos_offset_)
        // We need to remove pos_offset_ to get the raw position within our internal counter
        // Window offset = (raw_pos - (current_pos - w))
        u32x8 raw_pos = min_pos - pos_offset_;
        u32x8 window_start = pos_ - u32x8::splat(static_cast<uint32_t>(w_));
        return raw_pos - window_start;
    }

private:
    static u32x8 simd_min(u32x8 a, u32x8 b) {
        if constexpr (LEFT) {
            return a.min(b);
        } else {
            return a.max(b);
        }
    }

    void compute_suffix_minima(u32x8 last_elem) {
        // Compute suffix minima in reverse
        u32x8 suffix_min = ring_buf_[w_ - 1];
        for (size_t i = w_ - 1; i > 0; --i) {
            suffix_min = simd_min(suffix_min, ring_buf_[i - 1]);
            ring_buf_[i - 1] = suffix_min;
        }

        // Reset prefix minimum to last element
        prefix_min_ = last_elem;
    }

    void reset_positions() {
        // Prevent overflow by resetting positions
        pos_ = pos_ - delta_;
        prefix_min_ = prefix_min_ - delta_;
        pos_offset_ = pos_offset_ + delta_;

        for (size_t i = 0; i < w_; ++i) {
            ring_buf_[i] = ring_buf_[i] - delta_;
        }
    }

    size_t w_ = 0;
    size_t total_len_ = 0;

    u32x8 val_mask_;
    u32x8 pos_mask_;
    u32x8 max_pos_;
    u32x8 delta_;
    u32x8 pos_;
    u32x8 pos_offset_;
    u32x8 prefix_min_;
    u32x8 fw_init_;

    RingBufSimd ring_buf_;
};

/// Convenience type aliases
using SlidingMinSimdLeft = SlidingMinSimd<true>;
using SlidingMinSimdRight = SlidingMinSimd<false>;

/// Scalar sliding minimum for syncmer detection
/// Tracks both minimum value and its position within the window
class SlidingMinScalar {
public:
    SlidingMinScalar() = default;

    void init(size_t w) {
        w_ = w;
        idx_ = 0;
        count_ = 0;
        buffer_.resize(w);
    }

    /// Push a hash value and return the position of minimum within window
    /// Returns (min_offset, is_valid) where min_offset is 0 for first, w-1 for last
    std::pair<size_t, bool> push(uint32_t hash) {
        // Store (hash, position within window)
        buffer_[idx_] = {hash, count_ % w_};
        idx_ = (idx_ + 1) % w_;
        count_++;

        if (count_ < w_) {
            return {0, false};  // Window not full yet
        }

        // Find minimum and its offset
        uint32_t min_hash = UINT32_MAX;
        size_t min_offset = 0;

        for (size_t i = 0; i < w_; ++i) {
            size_t pos = (idx_ + i) % w_;  // oldest to newest
            if (buffer_[pos].first < min_hash) {
                min_hash = buffer_[pos].first;
                min_offset = i;  // 0 = oldest (first), w-1 = newest (last)
            }
        }

        return {min_offset, true};
    }

    /// Check if current window is a closed syncmer
    /// (minimum is at first or last position)
    bool is_closed_syncmer() const {
        if (count_ < w_) return false;

        uint32_t min_hash = UINT32_MAX;
        size_t min_offset = 0;

        for (size_t i = 0; i < w_; ++i) {
            size_t pos = (idx_ + i) % w_;
            if (buffer_[pos].first < min_hash) {
                min_hash = buffer_[pos].first;
                min_offset = i;
            }
        }

        return (min_offset == 0) || (min_offset == w_ - 1);
    }

private:
    size_t w_ = 0;
    size_t idx_ = 0;
    size_t count_ = 0;
    std::vector<std::pair<uint32_t, size_t>> buffer_;
};

} // namespace csyncmer_simd
