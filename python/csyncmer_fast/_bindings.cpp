#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <iostream>

#include "csyncmer_fast/utils.h"
#include "csyncmer_fast/iterator_syng.h"

namespace py = pybind11;

// Wrapper class for the iterator

class PySyncmerIterator
{
private:
    SyncmerIteratorS * iterator;
public:
    
    // constructor
    PySyncmerIterator(std::string& sequence, size_t k, size_t s){
        iterator = syncmer_generator_createS(sequence.data(), sequence.size(), k, s);
        if (!iterator){
            throw std::runtime_error("Failed to create syncmer iterator. Check sequence and k-mer/s-mer length.");
        }
    }
    
    // Destroyer
    ~PySyncmerIterator(){
        if (iterator){
            syncmer_generator_destroyS(iterator);
        }
    }
    
    //Class is non-copyable
    PySyncmerIterator(const PySyncmerIterator&) = delete;
    PySyncmerIterator operator=(const PySyncmerIterator&) = delete;

    // Python iterator
    PySyncmerIterator& iter() {return *this;}

    py::tuple get_all_syncmers(){
        std::vector<U64> hash_values;
        std::vector<U64> kmer_position;
        std::vector<U64> smer_position;

        Syncmer64 syncmer = {0,0,0};
        bool keep {true};
        while(keep){
            
            if (syncmer_iterator_nextS(iterator, &syncmer)){
                hash_values.push_back(syncmer.hash_value);
                kmer_position.push_back(syncmer.kmer_position);
                smer_position.push_back(syncmer.smer_position);
            }
            
        }
        size_t num_elements = hash_values.size();
        std::cout << "SEEN " << num_elements << " SYNCMERS." << std::endl;
        py::array_t<U64> hash_array(num_elements);
        py::array_t<U64> kpos_array(num_elements);
        py::array_t<U64> spos_array(num_elements);
        
        std::memcpy(hash_array.mutable_data(), hash_values.data(), num_elements * sizeof(U64));
        std::memcpy(kpos_array.mutable_data(), kmer_position.data(), num_elements * sizeof(U64));
        std::memcpy(spos_array.mutable_data(), smer_position.data(), num_elements * sizeof(U64));

        return py::make_tuple(hash_array, kpos_array, spos_array);
    }

    py::object next() {
        Syncmer64 syncmer = {0,0,0};
        if (syncmer_iterator_nextS(iterator, &syncmer)){
            // auto [high, low] = split_uint128(syncmer.hash_value);
            return py::make_tuple(syncmer.hash_value, syncmer.kmer_position, syncmer.smer_position);
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
        .def(py::init<std::string&, size_t, size_t>(),
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
        .def("__next__", &PySyncmerIterator::next)
        .def("get_all_syncmers", &PySyncmerIterator::get_all_syncmers, "A function that returns all syncmer as 3 numpy vectors.");
}

