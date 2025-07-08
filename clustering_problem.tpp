#include "clustering_problem.hpp"

template<typename S>
ClusteringProblem<S>::ClusteringProblem(
    const std::vector<S> &samples,
    const std::function<int64_t(Utuple<3,S>)> &tripleCostCB,
    const std::function<int64_t(Utuple<2,S>)> &pairCostCB,
    std::string logFile
) : samples(samples) {
    // init the log file
    logs = std::ofstream(logFile);
    // init sample count
    sampleCount = samples.size();
    // init sample mapping
    sampleMapping.resize(sampleCount, {});
    for (int64_t i = 0; i < sampleCount; i++) {
        sampleMapping[i].push_back(i);
    }
    // init relevant pairs
    int64_t negativeSamples = 0, positiveSamples = 0;
    relevantPairs.resize(sampleCount, {});
    for (int64_t i = 0; i < sampleCount; i++) {
        for (int64_t j = i + 1; j < sampleCount; j++) {
            Upair indexPair({i, j});
            Utuple<2,S> samplePair({samples[i], samples[j]});
            int64_t c = pairCostCB(samplePair);
            if (c) {
                pairCosts[indexPair] = c;
                relevantPairs[i].push_back(j);
                relevantPairs[j].push_back(i);
                if (c > 0) positiveSamples++;
                if (c < 0) negativeSamples++;
            }
        }
    }
    // init relevant triples (preserve j < k)
    relevantTriples.resize(sampleCount, {});
    for (int64_t i = 0; i < sampleCount; i++) {
        for (int64_t j = i + 1; j < sampleCount; j++) {
            for (int64_t k = j + 1; k < sampleCount; k++) {
                Utriple indexTriple({i, j, k});
                Utuple<3,S> sampleTriple({samples[i], samples[j], samples[k]});
                int64_t c = tripleCostCB(sampleTriple);
                if (c) {
                    tripleCosts[indexTriple] = c;
                    relevantTriples[i].push_back(std::make_pair(j, k)); // j < k!
                    relevantTriples[j].push_back(std::make_pair(i, k)); // i < k!
                    relevantTriples[k].push_back(std::make_pair(i, j)); // j < k!
                    if (c > 0) positiveSamples++;
                    if (c < 0) negativeSamples++;
                }
            }
        }
    }
    logs<< "------------------------------------\n";
    logs<< "Clustering problem has been inited with " << tripleCosts.size() << " relevant triples";
    logs<< ", (negative: " << negativeSamples << ", positive: " << positiveSamples << ")";
    logs<< std::endl;
}

template<typename S>
int64_t ClusteringProblem<S>::getCost(int64_t i, int64_t j) {
    Upair indexPair({i, j});
    if (pairCosts.count(indexPair)) return pairCosts[indexPair];
    return 0;
}

template<typename S>
int64_t ClusteringProblem<S>::getCost(int64_t i, int64_t j, int64_t k) {
    Utriple indexTriple({i, j, k});
    if (tripleCosts.count(indexTriple)) return tripleCosts[indexTriple];
    return 0;
}

template<typename S>
void ClusteringProblem<S>::solve() {
    resultingCost = 0;
    solve(std::vector<bool>(sampleCount, true));
    printResults();
}

template<typename S>
void ClusteringProblem<S>::solve(const std::vector<bool> &relevant) {
    int64_t subSampleCount = std::count(std::begin(relevant), std::end(relevant), true);
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
std::map<Upair, int64_t> ClusteringProblem<S>::getLabels() {
    return label;
}

template<typename S>
bool ClusteringProblem<S>::isSolvedCompletely() {
    for (int64_t i = 0; i < sampleCount; i++) {
        for (int64_t j = i + 1; j < sampleCount; j++) {
            if (!label[Upair({i, j})]) return false;
        }
    }
    return true;
}

template<typename S>
int64_t ClusteringProblem<S>::getSolutionCost() {
    return resultingCost;
}

template<typename S>
void ClusteringProblem<S>::printLabeling() {
    logs<< "Labeling: " << std::endl;
    logs<< "  ";
    for (auto s : samples) {
        logs<< s << " ";
    }
    logs<< std::endl;
    for (int64_t i = 0; i < samples.size(); i++) {
        logs<< samples[i] << " ";
        for (int64_t j = 0; j < samples.size(); j++) {
            if (i == j) {
                logs<< "- ";
                continue;
            }
            int64_t l = label[Upair({i, j})];
            if (l == 2) {
                logs<< 1;
            } else if (l == 1) {
                logs<< 0;
            } else {
                logs<< "x";
            }
            logs<< " ";
        }
        logs<< std::endl;
    }
    logs<< std::endl;
}

template<typename S>
void ClusteringProblem<S>::printClustering() {
    std::vector<std::vector<int64_t>> clustering;
    for (int64_t i = 0; i < sampleCount; i++) {
        if (sampleMapping[i].empty()) continue;
        clustering.push_back(sampleMapping[i]);
    }
    int64_t clusterCount = clustering.size();
    //
    logs<< "Cluster labeling: " << std::endl;
    logs<< "  ";
    for (int64_t i = 0; i < clusterCount; i++) {
        logs<< i << " ";
    }
    logs<< std::endl;
    for (int64_t ci = 0; ci < clusterCount; ci++) {
        logs<< ci << " ";
        int64_t i = clustering[ci][0];
        for (int64_t cj = 0; cj < clusterCount; cj++) {
            if (ci == cj) {
                logs<< "- ";
                continue;
            }
            int64_t j = clustering[cj][0];
            if (label[Upair({i, j})] == 1) { // cut clusters
                logs<< "0";
            } else { // label[i][j] == 0 unknown 
                logs<< "x";
            }
            logs<< " ";
            if (cj >= 10) logs<< " ";
        }
        logs<< std::endl;
    }
    logs<< std::endl;
    //
    logs<< "Unjoinable cluster triples: " << std::endl;
    for (int64_t ci = 0; ci < clusterCount; ci++) {
        int64_t i = clustering[ci][0];
        for (int64_t cj = ci + 1; cj < clusterCount; cj++) {
            int64_t j = clustering[cj][0];
            for (int64_t ck = cj + 1; ck < clusterCount; ck++) {
                int64_t k = clustering[ck][0];
                Utriple indexTriple({i, j, k});
                if (cutTriples.count(indexTriple)) {
                    logs<< "(" << ci << "," << cj << "," << ck << ") ";
                }
            }
        }
    }
    logs<< "\n" << std::endl;
    //
    logs<< "Clustering: " << std::endl;
    for (int64_t clusterInd = 0; clusterInd < clusterCount; clusterInd++) {
        logs<< "Cluster " << clusterInd << ": "; 
        std::vector<S> clusterSamples;
        for (auto originalI : clustering[clusterInd]) {
            clusterSamples.push_back(samples[originalI]);
        }
        std::sort(std::begin(clusterSamples), std::end(clusterSamples));
        for (auto &s : clusterSamples) {
            logs<< s;
        }
        logs<< std::endl;
    }
    logs<< std::endl;
}
  
template<typename S>
void ClusteringProblem<S>::printResults() {
    logs<< "---------------------------------------" << std::endl;
    logs<< "(0: cut; 1: joint; x: unknown)\n" << std::endl;
    printLabeling();
    printClustering();
    bool completeSolution = isSolvedCompletely();
    logs<< std::endl;
    logs<< "Problem solved: " << ((completeSolution) ? "completely" : "partially") << std::endl; 
    logs<< "Cost: " << resultingCost << std::endl;
    logs<< "---------------------------------------" << std::endl;
}

template<typename S>
bool ClusteringProblem<S>::applyIndependentSubproblemCut(const std::vector<bool> &relevant) {
    // 3.1
    logs<< "Trying independent subproblem cut (3.1)..." << std::endl;;
    std::vector<std::vector<int64_t>> partition;
    std::vector<bool> processed(sampleCount, false);
    for (int64_t i = 0; i < sampleCount; i++) {
        if (!relevant[i]) continue; // skip unrelevant samples 
        if (processed[i]) continue;
        // start the BFS from i
        std::vector<int64_t> chosen = {i};
        std::queue <int64_t> q;
        processed[i] = true;
        q.push(i);
        while (!q.empty()) {
            int64_t current = q.front();
            q.pop();
            for (auto j : relevantPairs[current]) {
                int64_t c = getCost(current, j);
                if (c >= 0) continue; 
                if (!processed[j]) {
                    processed[j] = true;
                    chosen.push_back(j);
                    q.push(j);
                }
            }
            for (auto [j, k] : relevantTriples[current]) {
                int64_t c = getCost(current, j, k);
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
    logs<< "* Applying the independent subproblem cut (3.1)" << std::endl;
    for (const auto &subproblemIndices : partition) {
        logs<< "Subproblem: ";
        for (int64_t i : subproblemIndices) {
            for (int64_t originalI : sampleMapping[i]) {
                logs<< samples[originalI];
            }
            logs<< " ";
        }
        logs<< std::endl;
    }
    logs<< std::endl;
    // create, solve the independent subproblems
    for (auto subproblemIndices : partition) {
        std::vector<bool> indexSubset(sampleCount, false);
        for (int64_t i : subproblemIndices) {
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
    for (int64_t i = 0; i < sampleCount; i++) {
        if (!relevant[i]) continue;
        if (!indexSubset[i]) continue; // i is in the index subset of the relevant subset
        for (int64_t j = 0; j < sampleCount; j++) {
            if (!relevant[j]) continue;
            if (indexSubset[j]) continue; // j is not in the subset but in the relevant subset
            for (int64_t originalI : sampleMapping[i]) {
                for (int64_t originalJ : sampleMapping[j]) {
                    label[Upair({originalI, originalJ})] = 1;
                }
            }
        }
    }
    // filter the relevant pairs and triples
    for (int64_t i = 0; i < sampleCount; i++) {
        if (!relevant[i]) continue;
        // filter the relevant pairs
        std::vector<int64_t> filteredPairs;
        for (int64_t j : relevantPairs[i]) {
            if (indexSubset[i] xor indexSubset[j]) {
                if (i < j) pairCosts.erase(Upair({i, j}));
            } else {
                filteredPairs.push_back(j);
            }
        }
        relevantPairs[i] = filteredPairs;
        // filter the relevant triples
        std::vector<std::pair<int64_t,int64_t>> filteredTriples;
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
int64_t ClusteringProblem<S>::solveMinCutForIndexSubset(
    const std::vector<bool> &indexSubset, 
    bool takeNegativeCosts, 
    bool takePositiveCosts, 
    bool globalMinCut, 
    int64_t source, 
    const std::vector<int64_t> &sinks
) {
    // apply 4.2
    int64_t vertices = std::count(std::begin(indexSubset), std::end(indexSubset), true);
    std::vector<int64_t> indexMapping(sampleCount);
    for (int64_t i = 0, sampleNode = 0; i < sampleCount; i++) {
        if (indexSubset[i]) {
            indexMapping[i] = sampleNode++;
        }
    }
    // compute the adjacency matrix by transforming triples (the costs in the matrix are not divided by 2 to avoid floating numbers)
    std::vector<std::vector<int64_t>> adjMatrix(vertices, std::vector<int64_t>(vertices, 0));
    for (int64_t i = 0; i < indexSubset.size(); i++) {
        if (!indexSubset[i]) continue;
        for (int64_t j : relevantPairs[i]) {
            if (i > j) continue; // consider only (i < j)
            if (!indexSubset[j]) continue;
            auto indexPair = Upair({i, j}); // sorted indices
            int64_t c = getCost(i, j);
            if (c < 0 && !takeNegativeCosts) continue;
            if (c > 0 && !takePositiveCosts) continue;
            int64_t i_node = indexMapping[indexPair[0]]; 
            int64_t j_node = indexMapping[indexPair[1]];
            adjMatrix[i_node][j_node] += 2*abs(c);
        }
        for (auto [j, k] : relevantTriples[i]) {
            if (i > j || j > k) continue; // consider only (i < j < k)
            if (!indexSubset[j] || !indexSubset[k]) continue;
            auto indexTriple = Utriple({i, j, k}); // sorted indices
            int64_t c = getCost(i, j, k);
            if (c < 0 && !takeNegativeCosts) continue;
            if (c > 0 && !takePositiveCosts) continue;
            int64_t i_node = indexMapping[indexTriple[0]]; 
            int64_t j_node = indexMapping[indexTriple[1]];
            int64_t k_node = indexMapping[indexTriple[2]];
            adjMatrix[i_node][j_node] += abs(c);
            adjMatrix[i_node][k_node] += abs(c);
            adjMatrix[j_node][k_node] += abs(c); 
        }
    }
    // create adjacency list from the adjacency matrix
    std::vector<std::tuple<int64_t,int64_t,int64_t>> edges;
    int64_t sumCosts = 0;
    for (int64_t i_node = 0; i_node < vertices; i_node++) {
        for (int64_t j_node = i_node + 1; j_node < vertices; j_node++) {
            int64_t c = adjMatrix[i_node][j_node];
            if (c) {
                edges.push_back(std::make_tuple(i_node, j_node, c));
                sumCosts += c;
            }
        }
    } 
    // solve the MinCut problem
    if (edges.empty()) return 0;
    int64_t minCut;
    if (globalMinCut) {
        minCut = MinCut::solveGlobalMinCut(edges);
    } else {
        int64_t superSink = vertices;
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
    logs<< "Join: ";
    for (int64_t k = 0; k < sampleCount; k++) {
        if (!indexSubset[k]) continue;
        for (int64_t originalK : sampleMapping[k]) {
            logs<< samples[originalK];
        }
        logs<< " ";
    }
    logs<< '\n' << std::endl;
    // aplying 5.1
    std::vector<int64_t> joinSamples;
    for (int64_t i = 0; i < sampleCount; i++) {
        if (indexSubset[i]) joinSamples.push_back(i);
    }
    // fix the labels
    for (int64_t indI = 0; indI < joinSamples.size(); indI++) {
        for (int64_t indJ = indI + 1; indJ < joinSamples.size(); indJ++) {
            int64_t i = joinSamples[indI], j = joinSamples[indJ];
            for (int64_t originalI : sampleMapping[i]) {
                for (int64_t originalJ : sampleMapping[j]) {
                    label[Upair({originalI, originalJ})] = 2;
                }
            }
        }
    }
    // update the sample mapping (modify the problem state)
    int64_t jointIndex = joinSamples[0];
    std::vector<int64_t> originaljoint;
    for (int64_t i : joinSamples) {
        for (auto originalI : sampleMapping[i]) {
            originaljoint.push_back(originalI);
        }
        sampleMapping[i].clear();
    }
    sampleMapping[jointIndex] = originaljoint;
    // compute the costs to the joint subset as well as the inner joining cost (modify the problem state)
    std::map<Upair, int64_t> subPairCosts; // additional variable in order not to overwrite the values before reading them
    std::map<Utriple, int64_t> subTripleCosts; // additional variable in order not tot overwrite the values before reading them
    for (int64_t i = 0; i < sampleCount; i++) {
        if (!relevant[i]) continue;
        for (int64_t j : relevantPairs[i]) {
            if (i > j) continue; // consider only one direction (another combination j > i also occures)
            Upair oldIndexPair({i, j});
            int64_t c = getCost(i, j);
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
            int64_t c = getCost(i, j, k);
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
    for (int64_t i = 0; i < sampleCount; i++) { // clear the relevant pairs and triples for the overwrite
        if (!relevant[i]) continue;
        relevantPairs[i].clear();
        relevantTriples[i].clear();
    }
    for (const auto& [indexPair, c] : pairCosts) { // overwrite the relevant pairs
        int64_t i = indexPair[0], j = indexPair[1];
        if (c) {
            relevantPairs[i].push_back(j);
            relevantPairs[j].push_back(i);
        }
    }
    for (const auto& [indexTriple, c] : tripleCosts) {
        int64_t i = indexTriple[0], j = indexTriple[1], k = indexTriple[2];
        if (c) {
            relevantTriples[i].push_back({j, k}); // j < k
            relevantTriples[j].push_back({i, k}); // i < k
            relevantTriples[k].push_back({i, j}); // i < j
        }
    }
    // solve the subproblem only for the relevant elements with one element instead of the joint samples
    std::vector<bool> newRelevant(sampleCount, false);
    for (int64_t i = 0; i < sampleCount; i++) {
        if (relevant[i] && ! indexSubset[i]) newRelevant[i] = true;
    }
    newRelevant[jointIndex] = true;
    solve(newRelevant);
    return;
}

template<typename S>
bool ClusteringProblem<S>::checkSubsetJoinForIndexSubset(const std::vector<bool> &indexSubset) {   
    // compute rhs
    int64_t lhsLowerBound = 0; // avoid MinCut computation for lhs>rhs (in particular the edge cases with lhs=0, rhs<0 because of 3.1 applied before)
    int64_t singleRhs = 0, doubleRhs = 0;
    for (int64_t i = 0; i < sampleCount; i++) {
        if (!indexSubset[i]) continue; // i is in R
        for (auto j : relevantPairs[i]) {
            int64_t c = getCost(i, j);
            if (indexSubset[j]) {
                if (i < j) lhsLowerBound += c;
            } else if (c < 0) {
                singleRhs += c;
            }
        }
        for (auto [j, k] : relevantTriples[i]) {
            int64_t c = getCost(i, j, k);
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
    int64_t rhs = singleRhs + doubleRhs/2;
    if (lhsLowerBound > rhs) return false;
    int64_t lhs = -solveMinCutForIndexSubset(indexSubset, true, false, true);
    return (lhs <= rhs);
}

template<typename S>
bool ClusteringProblem<S>::applySubsetJoin(const std::vector<bool> &relevant) {
    // subset join condition 3.11
    // heuristically construct and check the candidate sets R for possible joining'
    logs<< "Trying subset join (3.11)..." << std::endl;
    for (int64_t i = 0; i < sampleCount; i++) {
        if (!relevant[i]) continue;
        for (int64_t j = i + 1; j < sampleCount; j++) {
            if (!relevant[j]) continue;
            Upair indexPair({i, j});
            if (getCost(i, j) > 0) continue;
            std::vector<bool> indexSubset(sampleCount, false); // R, doesn't have to be a connected component
            indexSubset[i] = indexSubset[j] = true;
            std::vector<bool> joinIndexSubset;
            if (checkSubsetJoinForIndexSubset(indexSubset)) joinIndexSubset = indexSubset;
            std::set<int64_t> candidates;
            for (int64_t k = 0; k < sampleCount; k++) {
                if (!relevant[k]) continue;
                if (!indexSubset[k]) {
                    candidates.insert(k);
                }
            }
            std::function<int64_t(int64_t)> computeOffset = [&](int64_t i) {
                // positive offset if not mergeable
                if (indexSubset[i]) return int64_t(1);
                int64_t offset = 0;
                for (auto j : relevantPairs[i]) {
                    if (!indexSubset[j]) continue;
                    int64_t c = getCost(i, j);
                    if (c > 0) {
                        return int64_t(1);
                    } else {
                        offset += c;
                    }
                }
                for (auto [j, k] : relevantTriples[i]) {
                    if (!indexSubset[j] || !indexSubset[k]) continue; // 2 of 3 triple elements are already in R
                    int64_t c = getCost(i, j, k);
                    if (c > 0) {
                        return int64_t(1);
                    } else {
                        offset += c;
                    }
                }
                return offset;
            };
            while(!candidates.empty()) {
                std::vector<int64_t> badCandidates;
                int64_t bestCandidate = -1, bestOffset = 0;
                for (auto k : candidates) {
                    int64_t offset = computeOffset(k);
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
                logs<< "* Applying the subset join (3.11)" << std::endl;
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
    logs<< "Trying pair join (3.4)..." << std::endl;
    for (int64_t i = 0; i < sampleCount; i++) {
        if (!relevant[i]) continue;
        for (int64_t j = i + 1; j < sampleCount; j++) {
            if (!relevant[j]) continue;
            // compute lhs
            int64_t lhs = 0;
            int64_t c = getCost(i, j);
            if (c < 0) lhs += -2*c;
            for (auto [k1, k2] : relevantTriples[i]) {
                if (k1 != j && k2 != j) continue; // k1 or k2 is j
                int64_t c = getCost(i, k1, k2);
                if (c < 0) lhs += -c;
            }
            // compute rhs (assume that the underlying graph is connected after applying the 3.1)
            int64_t rhs = solveMinCutForIndexSubset(relevant, true, true, false, i, {j});
            if (lhs >= rhs) {
                std::vector<bool> indexSubset(sampleCount, false);
                indexSubset[i] = indexSubset[j] = true;
                logs<< "* Applying the pair join (3.4)" << std::endl;
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
    logs<< "Trying complex pair join (3.6)..." << std::endl;
    for (int64_t i = 0; i < sampleCount; i++) {
        if (!relevant[i]) continue;
        for (int64_t j = 0; j < sampleCount; j++) {
            if (j == i) continue;
            if (!relevant[j]) continue;
            for (int64_t k = 0; k < sampleCount; k++) {
                if (k == i || k == j) continue;
                if (!relevant[k]) continue; 
                // compute lhs1, lhs2, lhs3
                int64_t lhs1 = 0, lhs2 = 0, lhs3 = 0;
                {
                    int64_t c = getCost(i, j, k);
                    lhs3 += c;
                    if (c < 0) {
                        lhs1 += -c;
                        lhs2 += -c; 
                    }
                }
                {
                    int64_t c = getCost(i, j);
                    lhs3 += c;
                    if (c < 0) lhs1 += -2*c; 
                }
                {
                    int64_t c = getCost(j, k);
                    lhs3 += c;
                    if (c < 0) lhs2 += -2*c; 
                }
                {
                    int64_t c = getCost(i, k);
                    lhs3 += c;
                    if (c < 0) {
                        lhs1 += -2*c;
                        lhs2 += -2*c;
                    }
                }
                for (auto [p, q] : relevantTriples[i]) {
                    if (p == j || p == k || q == j || q == k) {
                        int64_t c = getCost(i, p, q);
                        if (c < 0) lhs1 += -c;
                    }
                }
                for (auto [p, q] : relevantTriples[k]) {
                    if (p == j || p == i || q == j || q == i) {
                        int64_t c = getCost(k, p, q);
                        if (c < 0) lhs2 += -c;
                    }
                }
                // compute rhs3
                int64_t rhs3 = 0;
                for (auto [p, q] : relevantTriples[i]) {
                    if (p == j || p == k || q == j || q == k) continue;
                    int64_t c = getCost(i, p, q);
                    rhs3 -= abs(c);
                }
                for (auto [p, q] : relevantTriples[j]) {
                    if (p == i || p == k || q == i || q == k) continue;
                    int64_t c = getCost(j, p, q);
                    rhs3 -= abs(c);
                }
                for (auto [p, q] : relevantTriples[k]) {
                    if (p == i || p == j || q == i || q == j) continue;
                    int64_t c = getCost(k, p, q);
                    rhs3 -= abs(c);
                }
                for (auto p : relevantPairs[i]) {
                    if (p == j || p == k) continue;
                    int64_t c = getCost(i, p);
                    rhs3 -= abs(c);
                }
                for (auto p : relevantPairs[j]) {
                    if (p == i || p == k) continue;
                    int64_t c = getCost(j, p);
                    rhs3 -= abs(c);
                }
                for (auto p : relevantPairs[k]) {
                    if (p == i || p == j) continue;
                    int64_t c = getCost(k, p);
                    rhs3 -= abs(c);
                }
                // check the 3d condition
                if (!(lhs3 <= rhs3)) continue; 
                // check the 1st and the 2d condition
                int64_t rhs1 = solveMinCutForIndexSubset(relevant, true, true, false, i, {j, k});
                int64_t rhs2 = solveMinCutForIndexSubset(relevant, true, true, false, k, {i, j});
                if (lhs1 >= rhs1 && lhs2 >= rhs2) {
                    std::vector<bool> indexSubset(sampleCount, false);
                    indexSubset[i] = indexSubset[k] = true; // join i and k
                    logs<< "* Applying the complex pair join (3.6)" << std::endl;
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
    logs<< "Trying explicit pair join (3.8)..." << std::endl;
    for (int64_t i = 0; i < sampleCount; i++) {
        if (!relevant[i]) continue;
        for (int64_t j = i + 1; j < sampleCount; j++) {
            if (!relevant[j]) continue;
            int64_t lhs = getCost(i, j);
            // compute rhs
            int64_t rhs = 0;
            for (int64_t k : relevantPairs[i]) {
                if (k == j) continue;
                int64_t c = getCost(i, k);
                if (c < 0) rhs += c;
            }
            for (int64_t k : relevantPairs[j]) {
                if (k == i) continue;
                int64_t c = getCost(j, k);
                if (c < 0) rhs += c;
            }
            for (auto [k1, k2] : relevantTriples[i]) {
                int64_t c = getCost(i, k1, k2);
                if (c < 0) rhs += c;
            }
            for (auto [k1, k2] : relevantTriples[j]) {
                if (k1 == i || k2 == i) continue; // the triples (i, j, *) have already been considered in the previos for-loop
                int64_t c = getCost(j, k1, k2);
                if (c < 0) rhs += c;
            }
            // check the condition
            if (lhs <= rhs) {
                std::vector<bool> indexSubset(sampleCount, false);
                indexSubset[i] = indexSubset[j] = true;
                logs<< "* Applying the explicit pair join (3.8)" << std::endl;
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
    logs<< "Trying explicit pair join via triple (3.9)..." << std::endl;
    for (int64_t i = 0; i < sampleCount; i++) {
        if (!relevant[i]) continue;
        for (int64_t j = 0; j < sampleCount; j++) {
            if (j == i) continue;
            if (!relevant[j]) continue;
            for (int64_t k = 0; k < sampleCount; k++) {
                if (k == i || k == j) continue;
                if (!relevant[k]) continue; 
                // compute costs for the triple
                int64_t cIJ = getCost(i, j);
                int64_t cIK = getCost(i, k);
                int64_t cJK = getCost(j, k);
                int64_t cIJK = getCost(i, j, k);
                // compute rhs
                int64_t singleR = 0, doubleR = 0;
                std::vector<int64_t> elementsR = {i, j, k};
                std::vector<bool> indexSubset(sampleCount, false);
                indexSubset[i] = indexSubset[j] = indexSubset[k] = true;
                for (int64_t i1 : elementsR) {
                    for (int64_t j1 : relevantPairs[i1]) {
                        if (indexSubset[j1]) continue;
                        int64_t c = getCost(i1, j1);
                        if (c < 0) singleR += c;
                    }
                    for (auto [j1, k1] : relevantTriples[i1]) {
                        if (indexSubset[j1] && indexSubset[k1]) continue;
                        int64_t c = getCost(i1, j1, k1);
                        if (c > 0) continue;
                        if (indexSubset[j1] xor indexSubset[k1]) {
                            doubleR += c;
                        } else { // only i1 is in R
                            singleR += c;
                        }
                    }
                }
                int64_t rhs = singleR + doubleR/2;
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
                    logs<< "* Applying the complex pair join (3.9)" << std::endl;
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
    logs<< "Trying triple join (3.5)..." << std::endl;
    // compute lhs base value
    int64_t lhsBase = 0;
    for (int64_t i = 0; i < sampleCount; i++) {
        if (!relevant[i]) continue;
        for (int64_t j : relevantPairs[i]) {
            if (i > j) continue;
            int64_t c = getCost(i, j);
            if (c > 0) lhsBase -= c;
        }
        for (auto [j, k] : relevantTriples[i]) {
            if (i > j) continue;
            int64_t c = getCost(i, j, k);
            if (c > 0) lhsBase -= c;
        }
    }
    // iterate over all unordered triples (i, j, k)
    for (int64_t i = 0; i < sampleCount; i++) {
        if (!relevant[i]) continue;
        for (int64_t j = 0; j < sampleCount; j++) {
            if (j == i) continue;
            if (!relevant[j]) continue;
            for (int64_t k = 0; k < sampleCount; k++) {
                if (k == i || k == j) continue;
                if (!relevant[k]) continue; 
                // compute lhs
                int64_t lhs = lhsBase;
                int64_t cIJ = getCost(i, j);
                if (cIJ < 0) lhs += -2*cIJ;
                int64_t cIK = getCost(i, k);
                if (cIK < 0) lhs += -2*cIK;
                int64_t cJK = getCost(j, k);
                if (cJK < 0) lhs += -cJK;
                int64_t cIJK = getCost(i, j, k);
                if (cIJK < 0) lhs += -2*cIJK;
                lhs += std::min(std::min(int64_t(0), cIJ), std::min(cIK, cJK)); // min for 4 cases
                if (lhs < 0) continue; // because rhs >= 0
                // compute rhs
                int64_t rhs = solveMinCutForIndexSubset(relevant, true, false, false, i, {j, k});
                if (lhs >= rhs) {
                    std::vector<bool> indexSubset(sampleCount, false);
                    indexSubset[i] = indexSubset[j] = indexSubset[k] = true; // join ijk
                    logs<< "* Applying the triple join (3.5)" << std::endl;
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
    logs<< "Trying pair cuts (3.2)..." << std::endl;
    for (int64_t i = 0; i < sampleCount; i++) {
        if (!relevant[i]) continue;
        for (int64_t j = i + 1; j < sampleCount; j++) {
            if (!relevant[j]) continue;
            if (label[Upair({sampleMapping[i][0], sampleMapping[j][0]})]) continue; // clusters are already cut 
            // check the condition
            int64_t lhs = 0;
            {
                int64_t c = getCost(i, j);
                if (c > 0) lhs += c;
            }
            int64_t rhs = solveMinCutForIndexSubset(relevant, true, false, false, i, {j});
            if (lhs >= rhs) {
                // cut all the pairs of the original samples
                for (int64_t originalI : sampleMapping[i]) {
                    for (int64_t originalJ : sampleMapping[j]) {
                        label[Upair({originalI, originalJ})] = 1; 
                    }
                } 
                logs<< "* Applying the pair cut (3.2)" << std::endl;
                logs<< "Cut pair: ";
                for (int64_t originalI : sampleMapping[i]) {
                    logs<< samples[originalI];
                }
                logs<< " ";
                for (int64_t originalJ : sampleMapping[j]) {
                    logs<< samples[originalJ];
                }
                logs<< '\n' << std::endl;
            }
        }
    }
}

template<typename S>
void ClusteringProblem<S>::applyTripleCuts(const std::vector<bool> &relevant) {
    // 3.3
    logs<< "Trying triple cuts (3.3)..." << std::endl;
    for (int64_t i = 0; i < sampleCount; i++) {
        if (!relevant[i]) continue;
        for (int64_t j = 0; j < sampleCount; j++) {
            if (j == i) continue;
            if (!relevant[j]) continue;
            if (label[Upair({sampleMapping[i][0], sampleMapping[j][0]})]) continue; // clusters are already cut
            for (int64_t k = 0; k < sampleCount; k++) {
                if (k == i || k == j) continue;
                if (!relevant[k]) continue; 
                if (cutTriples.count(Utriple({i, j, k}))) continue;
                if (label[Upair({sampleMapping[i][0], sampleMapping[k][0]})]) continue; // clusters are already cut
                if (label[Upair({sampleMapping[j][0], sampleMapping[k][0]})]) continue; // clusters are already cut
                int64_t lhs = 0;
                int64_t cIJK = getCost(i, j, k);
                int64_t cIJ = getCost(i, j);
                int64_t cIK = getCost(i, k);
                if (cIJK > 0) lhs += cIJK;
                if (cIJ > 0) lhs += cIJ;
                if (cIK > 0) lhs += cIK;
                int64_t rhs = solveMinCutForIndexSubset(relevant, true, false, false, i, {j, k});
                if (lhs >= rhs) {
                    cutTriples.insert(Utriple{i, j, k});
                    logs<< "* Applying the triple cut (3.3)" << std::endl;
                    logs<< "Cut triple: ";
                    for (int64_t originalI : sampleMapping[i]) {
                        logs<< samples[originalI];
                    }
                    logs<< " ";
                    for (int64_t originalJ : sampleMapping[j]) {
                        logs<< samples[originalJ];
                    }
                    logs<< " ";
                    for (int64_t originalK : sampleMapping[k]) {
                        logs<< samples[originalK];
                    }
                    logs<< '\n' << std::endl;
                }
            }
        }
    }
}