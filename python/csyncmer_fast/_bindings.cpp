#include <pybind11/pybind11.h>
#include <string>


#include "csyncmer_fast/iterator.h"
#include "csyncmer_fast/nthash_wrapper.h"
#include "csyncmer_fast/utils.h"

namespace py = pybind11;

std::pair<uint64_t, uint64_t> split_uint128(U128 value) {
    uint64_t low = static_cast<uint64_t>(value);
    uint64_t high = static_cast<uint64_t>(value >> 64);
    return std::make_pair(high, low);
}

// Wrapper class for the iterator

class PySyncmerIterator
{
private:
    SyncmerIterator * iterator;
public:
    
    // constructor
    PySyncmerIterator(const std::string& sequence, size_t k, size_t s){
        iterator = syncmer_generator_create(sequence.c_str(), k, s);
        if (!iterator){
            throw std::runtime_error("Failed to create syncmer iterator. Check sequence and k-mer/s-mer length.");
        }
    }
    
    // Destroyer
    ~PySyncmerIterator(){
        if (iterator){
            syncmer_generator_destroy(iterator);
        }
    }
    
    //Class is non-copyable
    PySyncmerIterator(const PySyncmerIterator&) = delete;
    PySyncmerIterator operator=(const PySyncmerIterator&) = delete;

    // Python iterator
    PySyncmerIterator& iter() {return *this;}

    py::object next() {
        Syncmer128 syncmer = {0,0,0};
        if (syncmer_iterator_next(iterator, &syncmer)){
            auto [high, low] = split_uint128(syncmer.hash_value);
            return py::make_tuple(high, low, syncmer.kmer_position, syncmer.smer_position);
        } else {
            throw py::stop_iteration();
        }
    }
};

// PYBIND MODULE DELCARATION

PYBIND11_MODULE(_bindings, m){
    m.doc() = "Fast syncmer generation using nthash on single thread.";

    // Binding of the iterator class
    py::class_<PySyncmerIterator>(m,"SyncmerIterator")
        .def(py::init<const std::string&, size_t, size_t>(),
            py::arg("sequence"), py::arg("k"), py::arg("s"),
            "Create a syncmer iterator\n\n"
            "Parameters\n"
            "sequence: str\n"
            "   DNA sequence (ACGT)\n"
            "k: int\n"
            "   k-mer size\n"
            "s: int\n"
            "   s-mer size\n (must be < k)")
        .def("__iter__", &PySyncmerIterator::iter)
        .def("__next__", &PySyncmerIterator::next);
}

