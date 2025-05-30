#ifndef CLUSTERING_PROBLEM_HPP
#define CLUSTERING_PROBLEM_HPP

#include <iostream>
#include <map>
#include <vector>
#include <array>
#include <queue>
#include <functional>
#include <stdexcept>
#include <algorithm>

#include "utuple.hpp"

using Upair = Utuple<2>;
using Utriple = Utuple<3>;


template <typename S = int> // domain S of samples
class ClusteringProblem {

    const std::vector<S> &samples;
    int sampleCount;
    
    std::vector<std::vector<int>> sampleMapping; 

    std::vector<std::vector<std::pair<int, int>>> relevantTriples;
    std::map<Utriple, int> tripleCosts;

    std::vector<std::vector<int>> relevantPairs;
    std::map<Upair, int> pairCosts;

    std::vector<std::vector<bool>> labelFixed, labelValue;
    int resultingCost;

    void solve(const std::vector<bool> &indexSubset);

public:

    explicit ClusteringProblem(
        const std::vector<S> &samples,
        const std::function<int(Utuple<3,S>)> &tripleCostCB,
        const std::function<int(Utuple<2,S>)> &pairCostCB = [](Utuple<2,S> p)->int{return 0;}
    );
    
    void solve();
};

#include "clustering_problem.tpp"

#endif