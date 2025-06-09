#ifndef CLUSTERING_INSTANCE_HPP
#define CLUSTERING_INSTANCE_HPP

#include <iostream>
#include <vector>
#include <functional>

#include "utuple.hpp"


template <typename S = int64_t>
struct ClusteringInstance {
    std::vector<S> samples;
    std::function<int64_t(Utuple<3,S>)> cost;
    std::function<int64_t(Utuple<2,S>)> pairCost;
    
    ClusteringInstance(
        std::vector<S> samples,
        std::function<int64_t(Utuple<3,S>)> tripleCostCB,
        std::function<int64_t(Utuple<2,S>)> pairCostCB= [](Utuple<2,S> p)->int64_t{return 0;}
    ) : samples(samples), cost(tripleCostCB), pairCost(pairCostCB) {}
};


#endif