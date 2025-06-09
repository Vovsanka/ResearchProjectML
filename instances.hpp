#ifndef INSTANCES_HPP
#define INSTANCES_HPP

#include <iostream>
#include <vector>
#include <functional>

#include "utuple.hpp"
#include "space.hpp"


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


extern const ClusteringInstance<char> SIMPLE_INSTANCE;
extern const ClusteringInstance<char> MULTICLUSTER_INSTANCE;
extern const ClusteringInstance<char> PYRAMID_INSTANCE1;
extern const ClusteringInstance<char> PYRAMID_INSTANCE2; 
extern const ClusteringInstance<char> PYRAMID_INSTANCE_UNSOLVABLE;

#endif