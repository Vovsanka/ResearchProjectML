#ifndef INSTANCES_HPP
#define INSTANCES_HPP

#include <iostream>
#include <vector>
#include <functional>

#include "utuple.hpp"
#include "space.hpp"


template <typename S>
struct ClusteringInstance {
    std::vector<S> samples;
    std::function<int64_t(Utuple<3,S>)> cost;
    std::function<int64_t(Utuple<2,S>)> pairCost;
    // TODO: save correct labels to compare with!!!

    ClusteringInstance(
        std::vector<S> samples,
        std::function<int64_t(Utuple<3,S>)> tripleCostCB,
        std::function<int64_t(Utuple<2,S>)> pairCostCB= [](Utuple<2,S> p)->int64_t{return 0;}
    ) : samples(samples), cost(tripleCostCB), pairCost(pairCostCB) {}
};

template <typename U>
std::function<int64_t(U)> doubleToIntCostWrapper(
    const std::function<double(U)> &costCB,
    int64_t multiplyBy = 10
) {
    return [costCB, multiplyBy](U utuple) -> int64_t {
        double c = costCB(utuple);
        return std::round(c*multiplyBy);
    };
}

extern const ClusteringInstance<char> SIMPLE_INSTANCE;
extern const ClusteringInstance<char> MULTICLUSTER_INSTANCE;
extern const ClusteringInstance<char> PYRAMID_INSTANCE1;
extern const ClusteringInstance<char> PYRAMID_INSTANCE2; 
extern const ClusteringInstance<char> PYRAMID_INSTANCE_UNSOLVABLE;
extern const ClusteringInstance<Space::Point> CUBIC_SPACE_INSTANCE;

#endif