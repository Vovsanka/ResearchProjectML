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
void ClusteringProblem<S>::solve(const std::vector<bool> &relevant) {
    int subSampleCount = std::count(std::begin(relevant), std::end(relevant), true);
    if (!subSampleCount) throw std::runtime_error("Cannot solve a clustering problem with no samples!");
    if (subSampleCount == 1) return; // trivial problem
    // apply partial optimality conditions (solve the subproblems if needed)
    if (applyIndependentSubproblemCut(relevant)) return;
    // if (applySubsetJoin(relevant))) return;
    // if (applyPairJoin(relevant))) return;
    // if (applyComplexPairJoin(relevant))) return;
    // if (applyExplicitPairJoin(relevant))) return;
    // if (applyExplicitPairJoinViaTriple(relevant))) return;
    // if (applyTripleJoin(relevant))) return;
    // applyPairCuts(relevant));
    // applyTripleCuts(relevant));
    return;
}

template<typename S>
void ClusteringProblem<S>::solve() {
    labelFixed.resize(sampleCount, std::vector<bool>(sampleCount, false));
    labelValue.resize(sampleCount, std::vector<bool>(sampleCount, false));
    resultingCost = 0;
    solve(std::vector<bool>(sampleCount, true));
}

template<typename S>
bool ClusteringProblem<S>::applyIndependentSubproblemCut(const std::vector<bool> &relevant) {
    std::vector<std::vector<int>> partition;
    std::vector<bool> processed(sampleCount, false);
    for (int i = 0; i < sampleCount; i++) {
        if (!relevant[i]) continue; // skip unrelevant samples 
        if (processed[i]) continue;
        // start the BFS from i
        std::vector<int> chosen = {i};
        std::queue <int> q;
        processed[i] = true;
        q.push(i);
        while (!q.empty()) {
            int current = q.front();
            q.pop();
            for (auto j : relevantPairs[current]) {
                Upair indexPair({current, j});
                int c = pairCosts[indexPair];
                if (c >= 0) continue; 
                if (!processed[j]) {
                    processed[j] = true;
                    chosen.push_back(j);
                    q.push(j);
                }
            }
            for (auto [j, k] : relevantTriples[current]) {
                Utriple indexTriple({current, j, k});
                int c = tripleCosts[indexTriple];
                if (c >= 0) continue; 
                if (!processed[j]) {
                    processed[j] = true;
                    chosen.push_back(j);
                    q.push(j);
                }
                if (!processed[k]) {
                    processed[k] = true;
                    chosen.push_back(k);
                    q.push(k);
                }
            }
        }
        partition.push_back(chosen);
    }
    if (partition.size() == 1) return false; // no cuts => no smaller subproblems
    std::cout << "Applying the independent subproblem cut (proposition 3.1)\nCut: ";
    for (auto subproblemIndices : partition) {
        for (int i : subproblemIndices) {
            for (int originalI : sampleMapping[i]) {
                std::cout << samples[originalI];
            }
            std::cout << " ";
        }
        std::cout << "| ";
    }
    std::cout << "\n" << std::endl;
    // create, solve the independent subproblems
    for (auto subproblemIndices : partition) {
        std::vector<bool> indexSubset(sampleCount, false);
        for (int i : subproblemIndices) {
            indexSubset[i] = true;
        } 
        cutIndexSubset(relevant, indexSubset);
        solve(indexSubset); // solve only for the index subset
    }
    return true;
}

template<typename S>
void ClusteringProblem<S>::cutIndexSubset(const std::vector<bool> &relevant, const std::vector<bool> &indexSubset) {
    // fix the labels
    for (int i = 0; i < sampleCount; i++) {
        if (!indexSubset[i]) continue; // i is in the index subset of the relevant subset
        for (int j = 0; j < sampleCount; j++) {
            if (!relevant[j] || indexSubset[j]) continue; // j is not in the subset but in the relevant subset
            for (int originalI : sampleMapping[i]) {
                for (int originalJ : sampleMapping[j]) {
                    labelFixed[originalI][originalJ] = labelFixed[originalJ][originalI] = true;
                    labelValue[originalI][originalJ] = labelValue[originalJ][originalI] = false;
                }
            }
        }
    }
    // filter the relevant pairs
    for (int i = 0; i < sampleCount; i++) {
        if (!relevant[i]) continue;
        std::vector<int> filtered;
        for (int j : relevantPairs[i]) {
            if (indexSubset[i] xor indexSubset[j]) {
                if (i < j) pairCosts.erase(Upair({i, j}));
            } else {
                filtered.push_back(j);
            }
        }
        relevantPairs[i] = filtered;
    }
    // filter the relevant triples
    for (int i = 0; i < sampleCount; i++) {
        if (!relevant[i]) continue;
        std::vector<std::pair<int,int>> filtered;
        for (auto [j, k] : relevantTriples[i]) {
            if (
                (indexSubset[i] xor indexSubset[j]) ||
                (indexSubset[i] xor indexSubset[k]) ||
                (indexSubset[j] xor indexSubset[k])) {
                if (i < j) tripleCosts.erase(Utriple({i, j, k}));
            } else {
                filtered.push_back({j, k});
            }
        }
        relevantTriples[i] = filtered;
    }
}