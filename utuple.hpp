#ifndef UTUPLE_HPP
#define UTUPLE_HPP

#include <initializer_list>


template <int N, typename S = int> // domain S
class Utuple { // unordered tuple
    std::array<S, N> t;
public:
    Utuple(std::initializer_list<S> elements);
    S operator[](int ind) const; // only getter, not setter
    bool operator<(const Utuple<N,S> &other) const; // comparator to be a key in std::map
};

#include "utuple.tpp"

#endif