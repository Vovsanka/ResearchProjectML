#include "utuple.hpp"

template <int64_t N, typename S>
Utuple<N,S>::Utuple(std::initializer_list<S> elements) {
    int64_t ind = 0;
    for (auto &elem : elements) {
        t[ind++] = elem;
    }
    if (ind != N) throw std::runtime_error("Unordered tuple has a wrong size! "); 
    std::sort(std::begin(t), std::end(t));
    for (int64_t i = 1; i < N; i++) {
        if (t[i] == t[i - 1]) throw std::runtime_error("Unordered tuple contains non-unique elements! ");
    }
}

template <int64_t N, typename S>
S Utuple<N,S>::operator[](int64_t ind) const { // only getter, not setter
    if (!(0 <= ind && ind < N)) throw std::runtime_error("Unordered tuple is accessed with a wrong index! ");
    return t[ind];
}

template <int64_t N, typename S>
bool Utuple<N,S>::operator<(const Utuple<N,S> &other) const { // comparator to be a key in std::map
    for (int64_t i = 0; i < N - 1; i++) {
        if (t[i] < other[i]) return true;
        if (t[i] > other[i]) return false;
    }
    return t[N - 1] < other[N - 1];
}