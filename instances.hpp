#ifndef INSTANCES_HPP
#define INSTANCES_HPP

#include <iostream>
#include <vector>
#include <array>
#include <map>
#include <functional>

#include <Eigen/Dense>

#include "utuple.hpp"
#include "space.hpp"


template <typename S>
struct ClusteringInstance {
    std::vector<S> unlabeledSamples;
    std::vector<int64_t> actualClustering;
    std::function<int64_t(Utuple<3,S>)> cost;
    std::function<int64_t(Utuple<2,S>)> pairCost;
    std::map<Utuple<2,int64_t>,int64_t> label;

    ClusteringInstance(
        std::vector<std::pair<S,int64_t>> labeledSamples,
        std::function<int64_t(Utuple<3,S>)> tripleCostCB,
        std::function<int64_t(Utuple<2,S>)> pairCostCB= [](Utuple<2,S> p)->int64_t{return 0;}
    ) : cost(tripleCostCB), pairCost(pairCostCB) {
        this->actualClustering = std::vector<int64_t>(labeledSamples.size());
        for (int64_t i = 0; i < labeledSamples.size(); i++) {
            this->unlabeledSamples.push_back(labeledSamples[i].first);
            this->actualClustering[i] = labeledSamples[i].second;
        }
    }

};


ClusteringInstance<Space::Point> generateSpaceInstance(
    int64_t planeCount,
    int64_t pointsPerPlane,
    double maxDistance,
    double maxNoise
);

std::function<int64_t(Utuple<3,Space::Point>)> createSpaceCostFunction(
    const std::vector<Space::Point> &points,
    double maxDistance,
    double maxNoise
);

extern const ClusteringInstance<char> SIMPLE_INSTANCE;
extern const ClusteringInstance<char> MULTICLUSTER_INSTANCE;
extern const ClusteringInstance<char> PYRAMID_INSTANCE1;
extern const ClusteringInstance<char> PYRAMID_INSTANCE2; 
extern const ClusteringInstance<char> PYRAMID_INSTANCE_UNSOLVABLE;

#endif