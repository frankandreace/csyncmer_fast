#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <vector>

#include "csyncmer_fast.h"

namespace py = pybind11;

// Forward-only iterator: yields positions (int)
class PySyncmerIterator
{
private:
    CsyncmerIterator64* iterator;
public:
    PySyncmerIterator(const std::string& sequence, size_t k, size_t s) {
        iterator = csyncmer_iterator_create_64(sequence.data(), sequence.size(), k, s);
        if (!iterator) {
            throw std::runtime_error(
                "Failed to create syncmer iterator. "
                "Check that len(sequence) >= k and 0 < s < k.");
        }
    }

    ~PySyncmerIterator() {
        if (iterator) {
            csyncmer_iterator_destroy_64(iterator);
        }
    }

    PySyncmerIterator(const PySyncmerIterator&) = delete;
    PySyncmerIterator& operator=(const PySyncmerIterator&) = delete;

    PySyncmerIterator& iter() { return *this; }

    py::int_ next() {
        size_t pos;
        if (csyncmer_iterator_next_64(iterator, &pos)) {
            return py::int_(pos);
        }
        throw py::stop_iteration();
    }

    py::array_t<uint64_t> get_all_positions() {
        std::vector<uint64_t> positions;
        size_t pos;
        while (csyncmer_iterator_next_64(iterator, &pos)) {
            positions.push_back(static_cast<uint64_t>(pos));
        }
        py::array_t<uint64_t> arr(positions.size());
        std::memcpy(arr.mutable_data(), positions.data(),
                     positions.size() * sizeof(uint64_t));
        return arr;
    }
};

// Canonical iterator: yields (position, strand) tuples
class PyCanonicalSyncmerIterator
{
private:
    CsyncmerIteratorCanonical64* iterator;
public:
    PyCanonicalSyncmerIterator(const std::string& sequence, size_t k, size_t s) {
        iterator = csyncmer_iterator_create_canonical_64(
            sequence.data(), sequence.size(), k, s);
        if (!iterator) {
            throw std::runtime_error(
                "Failed to create canonical syncmer iterator. "
                "Check that len(sequence) >= k and 0 < s < k.");
        }
    }

    ~PyCanonicalSyncmerIterator() {
        if (iterator) {
            csyncmer_iterator_destroy_canonical_64(iterator);
        }
    }

    PyCanonicalSyncmerIterator(const PyCanonicalSyncmerIterator&) = delete;
    PyCanonicalSyncmerIterator& operator=(const PyCanonicalSyncmerIterator&) = delete;

    PyCanonicalSyncmerIterator& iter() { return *this; }

    py::tuple next() {
        size_t pos;
        int strand;
        if (csyncmer_iterator_next_canonical_64(iterator, &pos, &strand)) {
            return py::make_tuple(pos, strand);
        }
        throw py::stop_iteration();
    }

    py::tuple get_all_positions() {
        std::vector<uint64_t> positions;
        std::vector<uint8_t> strands;
        size_t pos;
        int strand;
        while (csyncmer_iterator_next_canonical_64(iterator, &pos, &strand)) {
            positions.push_back(static_cast<uint64_t>(pos));
            strands.push_back(static_cast<uint8_t>(strand));
        }
        py::array_t<uint64_t> pos_arr(positions.size());
        py::array_t<uint8_t> strand_arr(strands.size());
        std::memcpy(pos_arr.mutable_data(), positions.data(),
                     positions.size() * sizeof(uint64_t));
        std::memcpy(strand_arr.mutable_data(), strands.data(),
                     strands.size() * sizeof(uint8_t));
        return py::make_tuple(pos_arr, strand_arr);
    }
};

// Module-level count functions wrapping SIMD
static size_t count_syncmers(const std::string& sequence, size_t k, size_t s) {
    return csyncmer_twostack_simd_32_count(sequence.data(), sequence.size(), k, s);
}

static size_t count_syncmers_canonical(const std::string& sequence, size_t k, size_t s) {
    return csyncmer_twostack_simd_32_canonical_count(sequence.data(), sequence.size(), k, s);
}

PYBIND11_MODULE(_bindings, m) {
    m.doc() = "Fast closed syncmer detection using ntHash rolling hash.";

    py::class_<PySyncmerIterator>(m, "SyncmerIterator")
        .def(py::init<const std::string&, size_t, size_t>(),
            py::arg("sequence"), py::arg("k"), py::arg("s"),
            "Create a forward-only syncmer iterator.\n\n"
            "Yields syncmer positions (0-based) in the sequence.\n\n"
            "Parameters\n"
            "----------\n"
            "sequence : str\n"
            "    DNA sequence (ACGT)\n"
            "k : int\n"
            "    k-mer size\n"
            "s : int\n"
            "    s-mer size (must be < k)")
        .def("__iter__", &PySyncmerIterator::iter)
        .def("__next__", &PySyncmerIterator::next)
        .def("get_all_positions", &PySyncmerIterator::get_all_positions,
             "Return all syncmer positions as a numpy uint64 array.");

    py::class_<PyCanonicalSyncmerIterator>(m, "CanonicalSyncmerIterator")
        .def(py::init<const std::string&, size_t, size_t>(),
            py::arg("sequence"), py::arg("k"), py::arg("s"),
            "Create a canonical (strand-independent) syncmer iterator.\n\n"
            "Yields (position, strand) tuples where strand is 0=forward, 1=RC.\n\n"
            "Parameters\n"
            "----------\n"
            "sequence : str\n"
            "    DNA sequence (ACGT)\n"
            "k : int\n"
            "    k-mer size\n"
            "s : int\n"
            "    s-mer size (must be < k)")
        .def("__iter__", &PyCanonicalSyncmerIterator::iter)
        .def("__next__", &PyCanonicalSyncmerIterator::next)
        .def("get_all_positions", &PyCanonicalSyncmerIterator::get_all_positions,
             "Return all syncmer positions and strands as (numpy uint64 array, numpy uint8 array).");

    m.def("count_syncmers", &count_syncmers,
          py::arg("sequence"), py::arg("k"), py::arg("s"),
          "Count syncmers using SIMD-accelerated TWOSTACK (approximate, 32-bit hash).\n\n"
          "Returns the number of closed syncmers in the sequence.");

    m.def("count_syncmers_canonical", &count_syncmers_canonical,
          py::arg("sequence"), py::arg("k"), py::arg("s"),
          "Count canonical syncmers using SIMD-accelerated TWOSTACK (approximate, 32-bit hash).\n\n"
          "Returns the number of strand-independent closed syncmers in the sequence.");
}
