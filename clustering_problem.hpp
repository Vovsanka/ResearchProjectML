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
#include "min_cut.hpp"

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

    std::map<Upair, int> label;
    int resultingCost;

    void solve(const std::vector<bool> &relevant);

    bool applyIndependentSubproblemCut(const std::vector<bool> &relevant);

    void createSolveCutSubproblem(const std::vector<bool> &relevant, const std::vector<bool> &indexSubset);

    int solveMinCutForIndexSubset(
        const std::vector<bool> &indexSubset,
        bool takeNegativeCosts, bool takePositiveCosts,
        bool globalMinCut, 
        int source = 0, 
        const std::vector<int> &sinks = std::vector<int>({0})
    );

    void createSolveJoinSubproblem(const std::vector<bool> &relevant, const std::vector<bool> &indexSubset);

    bool checkSubsetJoinForIndexSubset(const std::vector<bool> &indexSubset);

    bool applySubsetJoin(const std::vector<bool> &relevant);

    bool applyPairJoin(const std::vector<bool> &relevant);

    bool applyComplexPairJoin(const std::vector<bool> &relevant);

    bool applyExplicitPairJoin(const std::vector<bool> &relevant);

    bool applyExplicitPairJoinViaTriple(const std::vector<bool> &relevant);

public:

    explicit ClusteringProblem(
        const std::vector<S> &samples,
        const std::function<int(Utuple<3,S>)> &tripleCostCB,
        const std::function<int(Utuple<2,S>)> &pairCostCB = [](Utuple<2,S> p)->int{return 0;}
    );

    std::map<Upair, int> getLabels();

    bool isSolvedCompletely();

    int getCost();

    void printLabeling();

    void printClustering();

    void printResults();
    
    void solve();
};

#include "clustering_problem.tpp"

#endif