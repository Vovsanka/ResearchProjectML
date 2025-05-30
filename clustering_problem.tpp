#include "clustering_problem.hpp"

template<typename S>
ClusteringProblem<S>::ClusteringProblem(
    const std::vector<S> &samples,
    const std::function<int(Utuple<3,S>)> &tripleCostCB,
    const std::function<int(Utuple<2,S>)> &pairCostCB
) : samples(samples) {
    // init sample count
    sampleCount = samples.size();
    // init sample mapping
    sampleMapping.resize(sampleCount, {});
    for (int i = 0; i < sampleCount; i++) {
        sampleMapping[i].push_back(i);
    }
    // init relevant pairs
    relevantPairs.resize(sampleCount, {});
    for (int i = 0; i < sampleCount; i++) {
        for (int j = i + 1; j < sampleCount; j++) {
            Upair indexPair({i, j});
            Utuple<2,S> samplePair({samples[i], samples[j]});
            int c = pairCostCB(samplePair);
            if (c) {
                pairCosts[indexPair] = c;
                relevantPairs[i].push_back(j);
                relevantPairs[j].push_back(i);
            }
        }
    }
    // init relevant triples
    relevantTriples.resize(sampleCount, {});
    for (int i = 0; i < sampleCount; i++) {
        for (int j = i + 1; j < sampleCount; j++) {
            for (int k = j + 1; k < sampleCount; k++) {
                Utriple indexTriple({i, j, k});
                Utuple<3,S> sampleTriple({samples[i], samples[j], samples[k]});
                int c = tripleCostCB(sampleTriple);
                if (c) {
                    tripleCosts[indexTriple] = c;
                    relevantTriples[i].push_back(std::make_pair(j, k)); // j < k
                    relevantTriples[j].push_back(std::make_pair(i, k)); // i < k
                    relevantTriples[k].push_back(std::make_pair(i, j)); // j < k
                }
            }
        }
    }
}


template<typename S>
void ClusteringProblem<S>::solve(const std::vector<bool> &indexSubset) {
    int subSampleCount = std::count(std::begin(indexSubset), std::end(indexSubset), true);
    if (!subSampleCount) throw std::runtime_error("Cannot solve a clustering problem with no samples!");
    if (subSampleCount == 1) return; // trivial problem
    // apply partial optimality conditions (solve the subproblems if needed)
    // if (applyIndependentSubproblemCut(indexSubset))) return;
    // if (applySubsetJoin(indexSubset))) return;
    // if (applyPairJoin(indexSubset))) return;
    // if (applyComplexPairJoin(indexSubset))) return;
    // if (applyExplicitPairJoin(indexSubset))) return;
    // if (applyExplicitPairJoinViaTriple(indexSubset))) return;
    // if (applyTripleJoin(indexSubset))) return;
    // applyPairCuts(indexSubset));
    // applyTripleCuts(indexSubset));
    return;
}

template<typename S>
void ClusteringProblem<S>::solve() {
    labelFixed.resize(sampleCount, std::vector<bool>(sampleCount, false));
    labelValue.resize(sampleCount, std::vector<bool>(sampleCount, false));
    resultingCost = 0;
    solve(std::vector<bool>(sampleCount, true));
}

void solve() {
    
}