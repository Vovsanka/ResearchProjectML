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
    // init relevant triples (preserve j < k)
    relevantTriples.resize(sampleCount, {});
    for (int i = 0; i < sampleCount; i++) {
        for (int j = i + 1; j < sampleCount; j++) {
            for (int k = j + 1; k < sampleCount; k++) {
                Utriple indexTriple({i, j, k});
                Utuple<3,S> sampleTriple({samples[i], samples[j], samples[k]});
                int c = tripleCostCB(sampleTriple);
                if (c) {
                    tripleCosts[indexTriple] = c;
                    relevantTriples[i].push_back(std::make_pair(j, k)); // j < k!
                    relevantTriples[j].push_back(std::make_pair(i, k)); // i < k!
                    relevantTriples[k].push_back(std::make_pair(i, j)); // j < k!
                }
            }
        }
    }
}

template<typename S>
void ClusteringProblem<S>::solve() {
    resultingCost = 0;
    solve(std::vector<bool>(sampleCount, true));
    printResults();
}

template<typename S>
void ClusteringProblem<S>::solve(const std::vector<bool> &relevant) {
    int subSampleCount = std::count(std::begin(relevant), std::end(relevant), true);
    if (subSampleCount == 1) return; // trivial problem
    // apply partial optimality conditions (solve the subproblems if needed)
    if (applyIndependentSubproblemCut(relevant)) return;
    if (applySubsetJoin(relevant)) return;
    if (applyPairJoin(relevant)) return;
    if (applyComplexPairJoin(relevant)) return;
    if (applyExplicitPairJoin(relevant)) return;
    if (applyExplicitPairJoinViaTriple(relevant)) return;
    if (applyTripleJoin(relevant)) return;
    applyPairCuts(relevant);
    applyTripleCuts(relevant);
    return;
}

template<typename S>
std::map<Upair, int> ClusteringProblem<S>::getLabels() {
    return label;
}

template<typename S>
bool ClusteringProblem<S>::isSolvedCompletely() {
    for (int i = 0; i < sampleCount; i++) {
        for (int j = i + 1; j < sampleCount; j++) {
            if (!label[Upair({i, j})]) return false;
        }
    }
    return true;
}

template<typename S>
int ClusteringProblem<S>::getSolutionCost() {
    return resultingCost;
}

template<typename S>
void ClusteringProblem<S>::printLabeling() {
    std::cout << "Labeling: " << std::endl;
    std::cout << "  ";
    for (auto s : samples) {
        std::cout << s << " ";
    }
    std::cout << std::endl;
    for (int i = 0; i < samples.size(); i++) {
        std::cout << samples[i] << " ";
        for (int j = 0; j < samples.size(); j++) {
            if (i == j) {
                std::cout << "- ";
                continue;
            }
            int l = label[Upair({i, j})];
            if (l == 2) {
                std::cout << 1;
            } else if (l == 1) {
                std::cout << 0;
            } else {
                std::cout << "x";
            }
            std::cout << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

template<typename S>
void ClusteringProblem<S>::printClustering() {
    std::vector<std::vector<int>> clustering;
    for (int i = 0; i < sampleCount; i++) {
        if (sampleMapping[i].empty()) continue;
        clustering.push_back(sampleMapping[i]);
    }
    int clusterCount = clustering.size();
    for (int clusterInd = 0; clusterInd < clusterCount; clusterInd++) {
        std::cout << "Cluster " << clusterInd << ": "; 
        for (auto originalI : clustering[clusterInd]) {
            std::cout << samples[originalI];
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
    std::cout << "Clustering: " << std::endl;
    std::cout << "  ";
    for (int i = 0; i < clusterCount; i++) {
        std::cout << i << " ";
    }
    std::cout << std::endl;
    for (int ci = 0; ci < clusterCount; ci++) {
        std::cout << ci << " ";
        int i = clustering[ci][0];
        for (int cj = 0; cj < clusterCount; cj++) {
            if (ci == cj) {
                std::cout << "- ";
                continue;
            }
            int j = clustering[cj][0];
            if (label[Upair({i, j})] == 1) { // cut clusters
                std::cout << "0";
            } else { // label[i][j] == 0 unknown 
                std::cout << "x";
            }
            std::cout << " ";
        }
        std::cout << std::endl;
    }
}
  
template<typename S>
void ClusteringProblem<S>::printResults() {
    std::cout << "---------------------------------------" << std::endl;
    std::cout << "(0: cut; 1: joint; x: unknown)\n" << std::endl;
    printLabeling();
    printClustering();
    bool completeSolution = isSolvedCompletely();
    std::cout << std::endl;
    std::cout << "Problem solved: " << ((completeSolution) ? "completely" : "partially") << std::endl; 
    std::cout << "Cost: " << resultingCost << std::endl;
    std::cout << "---------------------------------------" << std::endl;
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
    std::cout << "* Applying the independent subproblem cut (3.1)" << std::endl;
    for (const auto &subproblemIndices : partition) {
        std::cout << "Subproblem: ";
        for (int i : subproblemIndices) {
            for (int originalI : sampleMapping[i]) {
                std::cout << samples[originalI];
            }
            std::cout << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
    // create, solve the independent subproblems
    for (auto subproblemIndices : partition) {
        std::vector<bool> indexSubset(sampleCount, false);
        for (int i : subproblemIndices) {
            indexSubset[i] = true;
        } 
        createSolveCutSubproblem(relevant, indexSubset);
    }
    return true;
}

template<typename S>
void ClusteringProblem<S>::createSolveCutSubproblem(
    const std::vector<bool> &relevant,
    const std::vector<bool> &indexSubset
) {
    // fix the labels
    for (int i = 0; i < sampleCount; i++) {
        if (!relevant[i]) continue;
        if (!indexSubset[i]) continue; // i is in the index subset of the relevant subset
        for (int j = 0; j < sampleCount; j++) {
            if (!relevant[j]) continue;
            if (indexSubset[j]) continue; // j is not in the subset but in the relevant subset
            for (int originalI : sampleMapping[i]) {
                for (int originalJ : sampleMapping[j]) {
                    label[Upair({originalI, originalJ})] = 1;
                }
            }
        }
    }
    // filter the relevant pairs and triples
    for (int i = 0; i < sampleCount; i++) {
        if (!relevant[i]) continue;
        // filter the relevant pairs
        std::vector<int> filteredPairs;
        for (int j : relevantPairs[i]) {
            if (indexSubset[i] xor indexSubset[j]) {
                if (i < j) pairCosts.erase(Upair({i, j}));
            } else {
                filteredPairs.push_back(j);
            }
        }
        relevantPairs[i] = filteredPairs;
        // filter the relevant triples
        std::vector<std::pair<int,int>> filteredTriples;
        for (auto [j, k] : relevantTriples[i]) {
            if (
                (indexSubset[i] xor indexSubset[j]) ||
                (indexSubset[i] xor indexSubset[k]) ||
                (indexSubset[j] xor indexSubset[k])) {
                if (i < j) tripleCosts.erase(Utriple({i, j, k})); // i < j < k
            } else {
                filteredTriples.push_back({j, k}); // j < k by induction
            }
        }
        relevantTriples[i] = filteredTriples;
    }
    // solve the subproblem only for the index subset (relevant elements in the subproblem)
    solve(indexSubset);
    return;
}

template<typename S>
int ClusteringProblem<S>::solveMinCutForIndexSubset(
    const std::vector<bool> &indexSubset, 
    bool takeNegativeCosts, 
    bool takePositiveCosts, 
    bool globalMinCut, 
    int source, 
    const std::vector<int> &sinks
) {
    // apply 4.2
    int vertices = std::count(std::begin(indexSubset), std::end(indexSubset), true);
    std::vector<int> indexMapping(sampleCount);
    for (int i = 0, sampleNode = 0; i < sampleCount; i++) {
        if (indexSubset[i]) {
            indexMapping[i] = sampleNode++;
        }
    }
    // compute the adjacency matrix by transforming triples (the costs in the matrix are not divided by 2 to avoid floating numbers)
    std::vector<std::vector<int>> adjMatrix(vertices, std::vector(vertices, 0));
    for (int i = 0; i < indexSubset.size(); i++) {
        if (!indexSubset[i]) continue;
        for (int j : relevantPairs[i]) {
            if (i > j) continue; // consider only (i < j)
            if (!indexSubset[j]) continue;
            auto indexPair = Upair({i, j}); // sorted indices
            int c = pairCosts[indexPair];
            if (c < 0 && !takeNegativeCosts) continue;
            if (c > 0 && !takePositiveCosts) continue;
            int i_node = indexMapping[indexPair[0]]; 
            int j_node = indexMapping[indexPair[1]];
            adjMatrix[i_node][j_node] += 2*abs(c);
        }
        for (auto [j, k] : relevantTriples[i]) {
            if (i > j || j > k) continue; // consider only (i < j < k)
            if (!indexSubset[j] || !indexSubset[k]) continue;
            auto indexTriple = Utriple({i, j, k}); // sorted indices
            int c = tripleCosts[indexTriple];
            if (c < 0 && !takeNegativeCosts) continue;
            if (c > 0 && !takePositiveCosts) continue;
            int i_node = indexMapping[indexTriple[0]]; 
            int j_node = indexMapping[indexTriple[1]];
            int k_node = indexMapping[indexTriple[2]];
            adjMatrix[i_node][j_node] += abs(c);
            adjMatrix[i_node][k_node] += abs(c);
            adjMatrix[j_node][k_node] += abs(c); 
        }
    }
    // create adjacency list from the adjacency matrix
    std::vector<std::tuple<int,int,int>> edges;
    int sumCosts = 0;
    for (int i_node = 0; i_node < vertices; i_node++) {
        for (int j_node = i_node + 1; j_node < vertices; j_node++) {
            int c = adjMatrix[i_node][j_node];
            if (c) {
                edges.push_back(std::make_tuple(i_node, j_node, c));
                sumCosts += c;
            }
        }
    } 
    // solve the MinCut problem
    int minCut;
    if (globalMinCut) {
        minCut = MinCut::solveGlobalMinCut(edges);
    } else {
        int superSink = vertices;
        for (auto sink : sinks) {
            edges.push_back(std::make_tuple(indexMapping[sink], superSink, sumCosts));
        }
        minCut = MinCut::solveMinCut(vertices + 1, edges, indexMapping[source], superSink);
    }
    // divide the MinCut result by 2 because all triples have been transformed with a double cost
    return minCut/2;
}

template<typename S>
void ClusteringProblem<S>::createSolveJoinSubproblem(const std::vector<bool> &relevant, const std::vector<bool> &indexSubset) {
    std::cout << "Join: ";
    for (int k = 0; k < sampleCount; k++) {
        if (!indexSubset[k]) continue;
        for (int originalK : sampleMapping[k]) {
            std::cout << samples[originalK];
        }
        std::cout << " ";
    }
    std::cout << '\n' << std::endl;
    // aplying 5.1
    std::vector<int> joinSamples;
    for (int i = 0; i < sampleCount; i++) {
        if (indexSubset[i]) joinSamples.push_back(i);
    }
    // fix the labels
    for (int indI = 0; indI < joinSamples.size(); indI++) {
        for (int indJ = indI + 1; indJ < joinSamples.size(); indJ++) {
            int i = joinSamples[indI], j = joinSamples[indJ];
            for (int originalI : sampleMapping[i]) {
                for (int originalJ : sampleMapping[j]) {
                    label[Upair({originalI, originalJ})] = 2;
                }
            }
        }
    }
    // update the sample mapping (modify the problem state)
    int jointIndex = joinSamples[0];
    std::vector<int> originalJoint;
    for (int i : joinSamples) {
        for (auto originalI : sampleMapping[i]) {
            originalJoint.push_back(originalI);
        }
        sampleMapping[i].clear();
    }
    sampleMapping[jointIndex] = originalJoint;
    // compute the costs to the joint subset as well as the inner joining cost (modify the problem state)
    std::map<Upair, int> subPairCosts; // additional variable in order not to overwrite the values before reading them
    std::map<Utriple, int> subTripleCosts; // additional variable in order not tot overwrite the values before reading them
    for (int i = 0; i < sampleCount; i++) {
        if (!relevant[i]) continue;
        for (int j : relevantPairs[i]) {
            if (i > j) continue; // consider only one direction (another combination j > i also occures)
            Upair oldIndexPair({i, j});
            int c = pairCosts[oldIndexPair];
            if (indexSubset[i] && indexSubset[j]) { // inner joining
                resultingCost += c;
                pairCosts.erase(oldIndexPair);
            } else if (indexSubset[i] || indexSubset[j]) { // modify the costs to the samples being joint
                if (indexSubset[i]) subPairCosts[Upair({j, jointIndex})] += c;
                if (indexSubset[j]) subPairCosts[Upair({i, jointIndex})] += c;
                pairCosts.erase(oldIndexPair);
            } // else: outer costs are left unchanged
        }
        for (auto [j, k] : relevantTriples[i]) {
            if (i > j) continue; // consider only one direction (another triples will also occur)
            Utriple oldIndexTriple({i, j, k});
            int c = tripleCosts[oldIndexTriple];
            if (indexSubset[i] && indexSubset[j] && indexSubset[k]) {
                resultingCost += c;
                tripleCosts.erase(oldIndexTriple);
            } else if (indexSubset[i] || indexSubset[j] || indexSubset[k]) { // modify the cost to the samples being joint
                if (indexSubset[i] && indexSubset[j]) subPairCosts[Upair({k, jointIndex})] += c;
                if (indexSubset[i] && indexSubset[k]) subPairCosts[Upair({j, jointIndex})] += c;
                if (indexSubset[j] && indexSubset[k]) subPairCosts[Upair({i, jointIndex})] += c;
                if (!indexSubset[i] && !indexSubset[j]) subTripleCosts[Utriple({i, j, jointIndex})] += c;
                if (!indexSubset[i] && !indexSubset[k]) subTripleCosts[Utriple({i, k, jointIndex})] += c;
                if (!indexSubset[j] && !indexSubset[k]) subTripleCosts[Utriple({j, k, jointIndex})] += c;
                tripleCosts.erase(oldIndexTriple);
            } // else: outer costs are left unchanged
        }
    }
    for (const auto& [indexPair, c] : subPairCosts) { // apply the modifications after all values have been read
        pairCosts[indexPair] = c; 
    }
    for (const auto& [indexTriple, c] : subTripleCosts) {// apply the modifications after all values have been read
        tripleCosts[indexTriple] = c;
    }
    for (int i = 0; i < sampleCount; i++) { // clear the relevant pairs and triples for the overwrite
        if (!relevant[i]) continue;
        relevantPairs[i].clear();
        relevantTriples[i].clear();
    }
    for (const auto& [indexPair, c] : pairCosts) { // overwrite the relevant pairs
        int i = indexPair[0], j = indexPair[1];
        if (c) {
            relevantPairs[i].push_back(j);
            relevantPairs[j].push_back(i);
        }
    }
    for (const auto& [indexTriple, c] : tripleCosts) {
        int i = indexTriple[0], j = indexTriple[1], k = indexTriple[2];
        if (c) {
            relevantTriples[i].push_back({j, k}); // j < k
            relevantTriples[j].push_back({i, k}); // i < k
            relevantTriples[k].push_back({i, j}); // i < j
        }
    }
    // solve the subproblem only for the relevant elements with one element instead of the joint samples
    std::vector<bool> newRelevant(sampleCount, false);
    for (int i = 0; i < sampleCount; i++) {
        if (relevant[i] && ! indexSubset[i]) newRelevant[i] = true;
    }
    newRelevant[jointIndex] = true;
    solve(newRelevant);
    return;
}

template<typename S>
bool ClusteringProblem<S>::checkSubsetJoinForIndexSubset(const std::vector<bool> &indexSubset) {   
    // compute rhs
    int lhsLowerBound = 0; // avoid MinCut computation for lhs>rhs (in particular the edge cases with lhs=0, rhs<0 because of 3.1 applied before)
    int singleRhs = 0, doubleRhs = 0;
    for (int i = 0; i < sampleCount; i++) {
        if (!indexSubset[i]) continue; // i is in R
        for (auto j : relevantPairs[i]) {
            int c = pairCosts[Upair({i, j})];
            if (indexSubset[j]) {
                if (i < j) lhsLowerBound += c;
            } else if (c < 0) {
                singleRhs += c;
            }
        }
        for (auto [j, k] : relevantTriples[i]) {
            int c = tripleCosts[Utriple({i, j, k})];
            if (indexSubset[j] && indexSubset[k]) {
                if (i < j) lhsLowerBound += c;
                continue; // j or k must be not in R
            }
            if (c > 0) continue; // omit positive costs
            if (!indexSubset[k] && !indexSubset[j]) {
                singleRhs += c;
            } else {
                doubleRhs += c; // since two elements of the triple are in R, the cost will be added twice
            }
        }
    }
    int rhs = singleRhs + doubleRhs/2;
    if (lhsLowerBound > rhs) return false;
    int lhs = -solveMinCutForIndexSubset(indexSubset, true, false, true);
    return (lhs <= rhs);
}

template<typename S>
bool ClusteringProblem<S>::applySubsetJoin(const std::vector<bool> &relevant) {
    // subset join condition 3.11
    // heuristically construct and check the candidate sets R for possible joining
    for (int i = 0; i < sampleCount; i++) {
        if (!relevant[i]) continue;
        for (int j = i + 1; j < sampleCount; j++) {
            if (!relevant[j]) continue;
            Upair indexPair({i, j});
            if (pairCosts.count(indexPair) && pairCosts[indexPair] > 0) continue;
            std::vector<bool> indexSubset(sampleCount, false); // R, doesn't have to be a connected component
            indexSubset[i] = indexSubset[j] = true;
            std::vector<bool> joinIndexSubset;
            if (checkSubsetJoinForIndexSubset(indexSubset)) joinIndexSubset = indexSubset;
            std::set<int> candidates;
            for (int k = 0; k < sampleCount; k++) {
                if (!relevant[k]) continue;
                if (!indexSubset[k]) {
                    candidates.insert(k);
                }
            }
            std::function<int(int)> computeOffset = [&](int i) {
                // positive offset if not mergeable
                if (indexSubset[i]) return 1;
                int offset = 0;
                for (auto j : relevantPairs[i]) {
                    if (!indexSubset[j]) continue;
                    int c = pairCosts[Upair({i, j})];
                    if (c > 0) {
                        return 1;
                    } else {
                        offset += c;
                    }
                }
                for (auto [j, k] : relevantTriples[i]) {
                    if (!indexSubset[j] || !indexSubset[k]) continue; // 2 of 3 triple elements are already in R
                    int c = tripleCosts[Utriple({i, j, k})];
                    if (c > 0) {
                        return 1;
                    } else {
                        offset += c;
                    }
                }
                return offset;
            };
            while(!candidates.empty()) {
                std::vector<int> badCandidates;
                int bestCandidate = -1, bestOffset = 0;
                for (auto k : candidates) {
                    int offset = computeOffset(k);
                    if (offset > 0) {
                        badCandidates.push_back(k);
                    }
                    if (offset <= bestOffset) {
                        bestOffset = offset;
                        bestCandidate = k;
                    }
                }
                for (auto k : badCandidates) {
                    candidates.erase(k);
                }
                if (bestCandidate != -1) {
                    indexSubset[bestCandidate] = true;
                    candidates.erase(bestCandidate);
                    if (checkSubsetJoinForIndexSubset(indexSubset)) joinIndexSubset = indexSubset;
                }
            }
            if (!joinIndexSubset.empty()) {
                std::cout << "* Applying the subset join (3.11)" << std::endl;
                createSolveJoinSubproblem(relevant, joinIndexSubset);
                return true;
            }
        }
    }
    return false;
}

template<typename S>
bool ClusteringProblem<S>::applyPairJoin(const std::vector<bool> &relevant) {
    // 3.4
    for (int i = 0; i < sampleCount; i++) {
        if (!relevant[i]) continue;
        for (int j = i + 1; j < sampleCount; j++) {
            if (!relevant[j]) continue;
            // compute lhs
            int lhs = 0;
            Upair indexPair({i, j});
            int c = pairCosts[indexPair];
            if (c < 0) lhs += -2*c;
            for (auto [k1, k2] : relevantTriples[i]) {
                if (k1 != j && k2 != j) continue; // k1 or k2 is j
                int c = tripleCosts[Utriple({i, k1, k2})];
                if (c < 0) lhs += -c;
            }
            // compute rhs (assume that the underlying graph is connected after applying the 3.1)
            int rhs = solveMinCutForIndexSubset(relevant, true, true, false, i, {j});
            if (lhs >= rhs) {
                std::vector<bool> indexSubset(sampleCount, false);
                indexSubset[i] = indexSubset[j] = true;
                std::cout << "* Applying the pair join (3.4)" << std::endl;
                createSolveJoinSubproblem(relevant, indexSubset);
                return true;
            }
        }
    }
    return false;
}

template<typename S>
bool ClusteringProblem<S>::applyComplexPairJoin(const std::vector<bool> &relevant) {
    // 3.6
    for (int i = 0; i < sampleCount; i++) {
        if (!relevant[i]) continue;
        for (int j = 0; j < sampleCount; j++) {
            if (j == i) continue;
            if (!relevant[j]) continue;
            for (int k = 0; k < sampleCount; k++) {
                if (k == i || k == j) continue;
                if (!relevant[k]) continue; 
                // compute lhs1, lhs2, lhs3
                int lhs1 = 0, lhs2 = 0, lhs3 = 0;
                Utriple indexTriple({i, j, k});
                {
                    int c = tripleCosts[indexTriple];
                    lhs3 += c;
                    if (c < 0) {
                        lhs1 += -c;
                        lhs2 += -c; 
                    }
                }
                Upair indexPairIJ({i, j}), indexPairIK({i, k}), indexPairJK({j, k});
                {
                    int c = pairCosts[indexPairIJ];
                    lhs3 += c;
                    if (c < 0) lhs1 += -2*c; 
                }
                {
                    int c = pairCosts[indexPairJK];
                    lhs3 += c;
                    if (c < 0) lhs2 += -2*c; 
                }
                {
                    int c = pairCosts[indexPairIK];
                    lhs3 += c;
                    if (c < 0) {
                        lhs1 += -2*c;
                        lhs2 += -2*c;
                    }
                }
                for (auto [p, q] : relevantTriples[i]) {
                    if (p == j || p == k || q == j || q == k) {
                        int c = tripleCosts[Utriple({i, p, q})];
                        if (c < 0) lhs1 += -c;
                    }
                }
                for (auto [p, q] : relevantTriples[k]) {
                    if (p == j || p == i || q == j || q == i) {
                        int c = tripleCosts[Utriple({k, p, q})];
                        if (c < 0) lhs2 += -c;
                    }
                }
                // compute rhs3
                int rhs3 = 0;
                for (auto [p, q] : relevantTriples[i]) {
                    if (p == j || p == k || q == j || q == k) continue;
                    int c = tripleCosts[Utriple({i, p, q})];
                    rhs3 -= abs(c);
                }
                for (auto [p, q] : relevantTriples[j]) {
                    if (p == i || p == k || q == i || q == k) continue;
                    int c = tripleCosts[Utriple({j, p, q})];
                    rhs3 -= abs(c);
                }
                for (auto [p, q] : relevantTriples[k]) {
                    if (p == i || p == j || q == i || q == j) continue;
                    int c = tripleCosts[Utriple({k, p, q})];
                    rhs3 -= abs(c);
                }
                for (auto p : relevantPairs[i]) {
                    if (p == j || p == k) continue;
                    int c = pairCosts[Upair({i, p})];
                    rhs3 -= abs(c);
                }
                for (auto p : relevantPairs[j]) {
                    if (p == i || p == k) continue;
                    int c = pairCosts[Upair({j, p})];
                    rhs3 -= abs(c);
                }
                for (auto p : relevantPairs[k]) {
                    if (p == i || p == j) continue;
                    int c = pairCosts[Upair({k, p})];
                    rhs3 -= abs(c);
                }
                // check the 3d condition
                if (!(lhs3 <= rhs3)) continue; 
                // check the 1st and the 2d condition
                int rhs1 = solveMinCutForIndexSubset(relevant, true, true, false, i, {j, k});
                int rhs2 = solveMinCutForIndexSubset(relevant, true, true, false, k, {i, j});
                if (lhs1 >= rhs1 && lhs2 >= rhs2) {
                    std::vector<bool> indexSubset(sampleCount, false);
                    indexSubset[i] = indexSubset[k] = true; // join i and k
                    std::cout << "* Applying the complex pair join (3.6)" << std::endl;
                    createSolveJoinSubproblem(relevant, indexSubset);
                    return true;
                }
            }
        }
    }
    return false;
}

template<typename S>
bool ClusteringProblem<S>::applyExplicitPairJoin(const std::vector<bool> &relevant) {
    // 3.8
    for (int i = 0; i < sampleCount; i++) {
        if (!relevant[i]) continue;
        for (int j = i + 1; j < sampleCount; j++) {
            if (!relevant[j]) continue;
            int lhs = 0;
            Upair indexPair({i, j});
            if (pairCosts[indexPair]) lhs += pairCosts[indexPair];
            // compute rhs
            int rhs = 0;
            for (int k : relevantPairs[i]) {
                if (k == j) continue;
                int c = pairCosts[Upair({i, k})];
                if (c < 0) rhs += c;
            }
            for (int k : relevantPairs[j]) {
                if (k == i) continue;
                int c = pairCosts[Upair({j, k})];
                if (c < 0) rhs += c;
            }
            for (auto [k1, k2] : relevantTriples[i]) {
                int c = tripleCosts[Utriple({i, k1, k2})];
                if (c < 0) rhs += c;
            }
            for (auto [k1, k2] : relevantTriples[j]) {
                if (k1 == i || k2 == i) continue; // the triples (i, j, *) have already been considered in the previos for-loop
                int c = tripleCosts[Utriple({j, k1, k2})];
                if (c < 0) rhs += c;
            }
            // check the condition
            if (lhs <= rhs) {
                std::vector<bool> indexSubset(sampleCount, false);
                indexSubset[i] = indexSubset[j] = true;
                std::cout << "* Applying the explicit pair join (3.8)" << std::endl;
                createSolveJoinSubproblem(relevant, indexSubset);
                return true;
            }
        }
    }
    return false;
}

template<typename S>
bool ClusteringProblem<S>::applyExplicitPairJoinViaTriple(const std::vector<bool> &relevant) {
    // 3.9
    for (int i = 0; i < sampleCount; i++) {
        if (!relevant[i]) continue;
        for (int j = 0; j < sampleCount; j++) {
            if (j == i) continue;
            if (!relevant[j]) continue;
            for (int k = 0; k < sampleCount; k++) {
                if (k == i || k == j) continue;
                if (!relevant[k]) continue; 
                // compute costs for the triple
                Upair indexPairIJ({i, j}), indexPairIK({i, k}), indexPairJK({j, k});
                Utriple indexTriple({i, j, k});
                int cIJ = pairCosts[indexPairIJ];
                int cIK = pairCosts[indexPairIK];
                int cJK = pairCosts[indexPairJK];
                int cIJK = tripleCosts[indexTriple];
                // compute rhs
                int singleR = 0, doubleR = 0;
                std::vector<int> elementsR = {i, j, k};
                std::vector<bool> indexSubset(sampleCount, false);
                indexSubset[i] = indexSubset[j] = indexSubset[k] = true;
                for (int i1 : elementsR) {
                    for (int j1 : relevantPairs[i1]) {
                        if (indexSubset[j1]) continue;
                        int c = pairCosts[Upair({i1, j1})];
                        if (c < 0) singleR += c;
                    }
                    for (auto [j1, k1] : relevantTriples[i1]) {
                        if (indexSubset[j1] && indexSubset[k1]) continue;
                        int c = tripleCosts[Utriple({i1, j1, k1})];
                        if (c > 0) continue;
                        if (indexSubset[j1] xor indexSubset[k1]) {
                            doubleR += c;
                        } else { // only i1 is in R
                            singleR += c;
                        }
                    }
                }
                int rhs = singleR + doubleR/2;
                // check the conditions
                if (
                    cIJ + cIK <= 0 &&
                    cIJ + cJK <= 0 &&
                    cIK + cJK <= 0 &&
                    cIJ + cIK + cJK <= 0 &&
                    2*cIJ + 2*cIK + 2*cJK + cIJK <= 0 &&
                    cIJ + cIK + cIJK <= rhs &&
                    cJK + cIK + cIJK <= rhs
                ) {
                    std::vector<bool> indexSubset(sampleCount, false);
                    indexSubset[i] = indexSubset[k] = true; // join i and k
                    std::cout << "* Applying the complex pair join (3.9)" << std::endl;
                    createSolveJoinSubproblem(relevant, indexSubset);
                    return true;
                }
            }
        }
    }
    return false;
}

template<typename S>
bool ClusteringProblem<S>::applyTripleJoin(const std::vector<bool> &relevant) {
    // 3.5
    // compute lhs base value
    int lhsBase = 0;
    for (int i = 0; i < sampleCount; i++) {
        if (!relevant[i]) continue;
        for (int j : relevantPairs[i]) {
            if (i > j) continue;
            int c = pairCosts[Upair({i, j})];
            if (c > 0) lhsBase -= c;
        }
        for (auto [j, k] : relevantTriples[i]) {
            if (i > j) continue;
            int c = tripleCosts[Utriple({i, j, k})];
            if (c > 0) lhsBase -= c;
        }
    }
    // iterate over all unordered triples (i, j, k)
    for (int i = 0; i < sampleCount; i++) {
        if (!relevant[i]) continue;
        for (int j = 0; j < sampleCount; j++) {
            if (j == i) continue;
            if (!relevant[j]) continue;
            for (int k = 0; k < sampleCount; k++) {
                if (k == i || k == j) continue;
                if (!relevant[k]) continue; 
                // compute lhs
                int lhs = lhsBase;
                Upair indexPairIJ({i, j}), indexPairIK({i, k}), indexPairJK({j, k});
                Utriple indexTriple({i, j, k});
                int cIJ = pairCosts[indexPairIJ];
                if (cIJ < 0) lhs += -2*cIJ;
                int cIK = pairCosts[indexPairIK];
                if (cIK < 0) lhs += -2*cIK;
                int cJK = pairCosts[indexPairJK];
                if (cJK < 0) lhs += -cJK;
                int cIJK = tripleCosts[indexTriple];
                if (cIJK < 0) lhs += -2*cIJK;
                lhs += std::min(std::min(0, cIJ), std::min(cIK, cJK)); // min for 4 cases
                // compute rhs
                int rhs = solveMinCutForIndexSubset(relevant, true, false, false, i, {j, k});
                if (lhs >= rhs) {
                    std::vector<bool> indexSubset(sampleCount, false);
                    indexSubset[i] = indexSubset[j] = indexSubset[k] = true; // join ijk
                    std::cout << "* Applying the triple join (3.5)" << std::endl;
                    createSolveJoinSubproblem(relevant, indexSubset);
                    return true;
                }
            }
        }
    }
    return false;
}

template<typename S>
void ClusteringProblem<S>::applyPairCuts(const std::vector<bool> &relevant) {
    // 3.2
    for (int i = 0; i < sampleCount; i++) {
        if (!relevant[i]) continue;
        for (int j = i + 1; j < sampleCount; j++) {
            if (!relevant[j]) continue;
            // check the condition
            int lhs = 0;
            Upair indexPair({i, j});
            {
                int c = pairCosts[indexPair];
                if (c > 0) lhs += c;
            }
            int rhs = solveMinCutForIndexSubset(relevant, true, false, false, i, {j});
            if (lhs >= rhs) {
                // cut all the pairs of the original samples
                for (int originalI : sampleMapping[i]) {
                    for (int originalJ : sampleMapping[j]) {
                        label[Upair({originalI, originalJ})] = 1; 
                    }
                } 
                std::cout << "* Applying the pair cut (3.2)" << std::endl;
                std::cout << "Cut pair: ";
                for (int originalI : sampleMapping[i]) {
                    std::cout << samples[originalI];
                }
                std::cout << " ";
                for (int originalJ : sampleMapping[j]) {
                    std::cout << samples[originalJ];
                }
                std::cout << '\n' << std::endl;
            }
        }
    }
}

template<typename S>
void ClusteringProblem<S>::applyTripleCuts(const std::vector<bool> &relevant) {
    for (int i = 0; i < sampleCount; i++) {
        if (!relevant[i]) continue;
        for (int j = 0; j < sampleCount; j++) {
            if (j == i) continue;
            if (!relevant[j]) continue;
            for (int k = 0; k < sampleCount; k++) {
                if (k == i || k == j) continue;
                if (!relevant[k]) continue; 
                if (cutTriples.count(Utriple({i, j, k}))) continue;
                int lhs = 0;
                int cIJK = tripleCosts[Utriple({i, j, k})];
                int cIJ = pairCosts[Upair({i, j})];
                int cIK = pairCosts[Upair({i, k})];
                if (cIJK > 0) lhs += cIJK;
                if (cIJ > 0) lhs += cIJ;
                if (cIK > 0) lhs += cIK;
                int rhs = solveMinCutForIndexSubset(relevant, true, false, false, i, {j, k});
                if (lhs >= rhs) {
                    cutTriples.insert(Utriple{i, j, k});
                    std::cout << "* Applying the triple cut (3.3)" << std::endl;
                    std::cout << "Cut triple: ";
                    for (int originalI : sampleMapping[i]) {
                        std::cout << samples[originalI];
                    }
                    std::cout << " ";
                    for (int originalJ : sampleMapping[j]) {
                        std::cout << samples[originalJ];
                    }
                    std::cout << " ";
                    for (int originalK : sampleMapping[k]) {
                        std::cout << samples[originalK];
                    }
                    std::cout << '\n' << std::endl;
                }
            }
        }
    }
}