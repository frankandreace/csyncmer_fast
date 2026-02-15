import numpy as np
import csyncmer_fast
from csyncmer_fast import (
    SyncmerIterator,
    CanonicalSyncmerIterator,
    count_syncmers,
    count_syncmers_canonical,
)

# Sequence long enough for k=15, s=8 (need at least 15 bases)
SEQUENCE = "ACGTACGTACGTACGTACGTACGTACGTACGT"  # 32 bases
K = 15
S = 8


class TestImport:
    def test_version(self):
        assert hasattr(csyncmer_fast, '__version__')

    def test_exports(self):
        assert hasattr(csyncmer_fast, 'SyncmerIterator')
        assert hasattr(csyncmer_fast, 'CanonicalSyncmerIterator')
        assert hasattr(csyncmer_fast, 'count_syncmers')
        assert hasattr(csyncmer_fast, 'count_syncmers_canonical')


class TestSyncmerIterator:
    def test_yields_integers(self):
        positions = list(SyncmerIterator(SEQUENCE, K, S))
        assert len(positions) > 0
        for pos in positions:
            assert isinstance(pos, int)

    def test_positions_in_range(self):
        num_kmers = len(SEQUENCE) - K + 1
        for pos in SyncmerIterator(SEQUENCE, K, S):
            assert 0 <= pos < num_kmers

    def test_get_all_positions(self):
        arr = SyncmerIterator(SEQUENCE, K, S).get_all_positions()
        assert isinstance(arr, np.ndarray)
        assert arr.dtype == np.uint64
        assert len(arr) > 0

    def test_iterator_matches_batch(self):
        positions_iter = list(SyncmerIterator(SEQUENCE, K, S))
        positions_batch = SyncmerIterator(SEQUENCE, K, S).get_all_positions()
        assert len(positions_iter) == len(positions_batch)
        for a, b in zip(positions_iter, positions_batch):
            assert a == b

    def test_invalid_params(self):
        # s >= k
        try:
            SyncmerIterator(SEQUENCE, 5, 5)
            assert False, "Should have raised"
        except RuntimeError:
            pass
        # sequence too short
        try:
            SyncmerIterator("ACGT", 15, 8)
            assert False, "Should have raised"
        except RuntimeError:
            pass

    def test_deterministic(self):
        a = list(SyncmerIterator(SEQUENCE, K, S))
        b = list(SyncmerIterator(SEQUENCE, K, S))
        assert a == b


class TestCanonicalSyncmerIterator:
    def test_yields_tuples(self):
        results = list(CanonicalSyncmerIterator(SEQUENCE, K, S))
        assert len(results) > 0
        for pos, strand in results:
            assert isinstance(pos, int)
            assert strand in (0, 1)

    def test_positions_in_range(self):
        num_kmers = len(SEQUENCE) - K + 1
        for pos, strand in CanonicalSyncmerIterator(SEQUENCE, K, S):
            assert 0 <= pos < num_kmers

    def test_get_all_positions(self):
        positions, strands = CanonicalSyncmerIterator(SEQUENCE, K, S).get_all_positions()
        assert isinstance(positions, np.ndarray)
        assert isinstance(strands, np.ndarray)
        assert positions.dtype == np.uint64
        assert strands.dtype == np.uint8
        assert len(positions) == len(strands)
        assert len(positions) > 0

    def test_iterator_matches_batch(self):
        results_iter = list(CanonicalSyncmerIterator(SEQUENCE, K, S))
        positions_batch, strands_batch = CanonicalSyncmerIterator(SEQUENCE, K, S).get_all_positions()
        assert len(results_iter) == len(positions_batch)
        for (pos_i, strand_i), pos_b, strand_b in zip(
            results_iter, positions_batch, strands_batch
        ):
            assert pos_i == pos_b
            assert strand_i == strand_b

    def test_deterministic(self):
        a = list(CanonicalSyncmerIterator(SEQUENCE, K, S))
        b = list(CanonicalSyncmerIterator(SEQUENCE, K, S))
        assert a == b


class TestCountFunctions:
    def test_count_syncmers(self):
        count = count_syncmers(SEQUENCE, K, S)
        assert isinstance(count, int)
        assert count > 0

    def test_count_canonical(self):
        count = count_syncmers_canonical(SEQUENCE, K, S)
        assert isinstance(count, int)
        assert count > 0

    def test_count_reasonable(self):
        """SIMD count uses 32-bit hash, iterator uses 64-bit.
        Different hash sizes produce different syncmer counts due to tie-breaking,
        especially on short/repetitive sequences."""
        iter_count = len(list(SyncmerIterator(SEQUENCE, K, S)))
        simd_count = count_syncmers(SEQUENCE, K, S)
        assert abs(iter_count - simd_count) <= max(2, iter_count // 3)

    def test_count_canonical_reasonable(self):
        iter_count = len(list(CanonicalSyncmerIterator(SEQUENCE, K, S)))
        simd_count = count_syncmers_canonical(SEQUENCE, K, S)
        assert abs(iter_count - simd_count) <= max(2, iter_count // 3)
