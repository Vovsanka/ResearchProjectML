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


template <typename S = int64_t> // domain S of samples
class ClusteringProblem {

    const std::vector<S> &samples;
    int64_t sampleCount;
    
    std::vector<std::vector<int64_t>> sampleMapping; 

    std::vector<std::vector<std::pair<int64_t, int64_t>>> relevantTriples;
    std::map<Utriple, int64_t> tripleCosts;

    std::vector<std::vector<int64_t>> relevantPairs;
    std::map<Upair, int64_t> pairCosts;

    std::map<Upair, int64_t> label;
    std::set<Utriple> cutTriples;
    int64_t resultingCost;

    void solve(const std::vector<bool> &relevant);

    bool applyIndependentSubproblemCut(const std::vector<bool> &relevant);

    void createSolveCutSubproblem(const std::vector<bool> &relevant, const std::vector<bool> &indexSubset);

    int64_t solveMinCutForIndexSubset(
        const std::vector<bool> &indexSubset,
        bool takeNegativeCosts, bool takePositiveCosts,
        bool globalMinCut, 
        int64_t source = 0, 
        const std::vector<int64_t> &sinks = std::vector<int64_t>({0})
    );

    int64_t getCost(int64_t i, int64_t j);

    int64_t getCost(int64_t i, int64_t j, int64_t k);

    void createSolveJoinSubproblem(const std::vector<bool> &relevant, const std::vector<bool> &indexSubset);

    bool checkSubsetJoinForIndexSubset(const std::vector<bool> &indexSubset);

    bool applySubsetJoin(const std::vector<bool> &relevant);

    bool applyPairJoin(const std::vector<bool> &relevant);

    bool applyComplexPairJoin(const std::vector<bool> &relevant);

    bool applyExplicitPairJoin(const std::vector<bool> &relevant);

    bool applyExplicitPairJoinViaTriple(const std::vector<bool> &relevant);

    bool applyTripleJoin(const std::vector<bool> &relevant);

    void applyPairCuts(const std::vector<bool> &relevant);

    void applyTripleCuts(const std::vector<bool> &relevant);

public:

    explicit ClusteringProblem(
        const std::vector<S> &samples,
        const std::function<int64_t(Utuple<3,S>)> &tripleCostCB,
        const std::function<int64_t(Utuple<2,S>)> &pairCostCB = [](Utuple<2,S> p)->int64_t{return 0;}
    );

    std::map<Upair, int64_t> getLabels();

    bool isSolvedCompletely();

    int64_t getSolutionCost();

    void print64_tLabeling();

    void print64_tClustering();

    void print64_tResults();
    
    void solve();
};

#include "clustering_problem.tpp"

#endif