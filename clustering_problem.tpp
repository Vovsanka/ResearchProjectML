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
int ClusteringProblem<S>::getResultingCost() {
    return resultingCost;
}

template<typename S>
void ClusteringProblem<S>::printResultingLabeling() {
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
}

template<typename S>
void ClusteringProblem<S>::printResultingClustering() {
    std::vector<std::vector<int>> cluster(sampleCount, std::vector<int>());
    std::vector<std::vector<int>> unknown(sampleCount, std::vector<int>());
    for (int i = 0; i < sampleCount; i++) {
        cluster[i].push_back(i);
        for (int j = i + 1; j < sampleCount; j++) {
            int l = label[Upair({i, j})];
            if (l == 2) {
                cluster[i].push_back(j);
                cluster[j].push_back(i);
            } else if (!l) {
                unknown[i].push_back(j);
                unknown[j].push_back(i);
            } 
        }
    }
    std::cout << "Clusters: " << std::endl; 
    std::vector<bool> shown(sampleCount, false);
    for (int i = 0; i < sampleCount; i++) {
        if (shown[i]) continue;
        for (auto j : cluster[i]) {
            std::cout << samples[j];
            shown[j] = true;
        }
        if (!unknown[i].empty()) std::cout << " ? ";
        for (auto j : unknown[i]) {
            std::cout << samples[j];
        }
        std::cout << std::endl;
    }
}
  
template<typename S>
void ClusteringProblem<S>::printResults() {
    printResultingLabeling();
    printResultingClustering();
    bool completeSolution = isSolvedCompletely();
    std::cout << std::endl;
    std::cout << "Problem solved: " << ((completeSolution) ? "completely" : "partially") << std::endl; 
    std::cout << "Resulting cost: " << resultingCost << std::endl;
}

template<typename S>
void ClusteringProblem<S>::solve() {
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
    std::cout << "Applying the independent subproblem cut (proposition 3.1)" << std::endl;
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
void ClusteringProblem<S>::cutIndexSubset(
    const std::vector<bool> &relevant,
    const std::vector<bool> &indexSubset
) {
    // fix the labels
    for (int i = 0; i < sampleCount; i++) {
        if (!relevant[i] || !indexSubset[i]) continue; // i is in the index subset of the relevant subset
        for (int j = 0; j < sampleCount; j++) {
            if (!relevant[j] || indexSubset[j]) continue; // j is not in the subset but in the relevant subset
            for (int originalI : sampleMapping[i]) {
                for (int originalJ : sampleMapping[j]) {
                    label[Upair({originalI, originalJ})] = label[Upair({originalJ, originalI})] = 1;
                }
            }
        }
    }
    // filter the relevant pairs
    for (int i = 0; i < sampleCount; i++) {
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
    return;
}

template<typename S>
int ClusteringProblem<S>::solveMinCutForIndexSubset(
    std::vector<bool> indexSubset, 
    bool takeNegativeCosts, 
    bool takePositiveCosts, 
    bool globalMinCut, 
    int source, 
    std::vector<int> sinks
) {
    // apply proposition 4.2
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