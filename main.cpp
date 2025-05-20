#include <iostream>
#include <set>
#include <map>
#include <vector>
#include <array>
#include <queue>
#include <functional>
#include <stdexcept>
#include <algorithm>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/push_relabel_max_flow.hpp>
#include <boost/graph/stoer_wagner_min_cut.hpp>
#include <boost/property_map/property_map.hpp>


std::pair<int, std::vector<bool>> solveMinCut(int vertices, std::vector<std::tuple<int,int,int>> edges, int s, int t) {
    // Bigraph!
    using namespace boost;

    // Define the graph type
    typedef adjacency_list<vecS, vecS, directedS,
        property<vertex_name_t, std::string>,
        property<edge_capacity_t, int,
            property<edge_residual_capacity_t, int,
                property<edge_reverse_t, adjacency_list<>::edge_descriptor>>>> Graph;
            
    Graph g(vertices); 
    auto capacity = get(edge_capacity, g);
    auto rev = get(edge_reverse, g);
    auto residual_capacity = get(edge_residual_capacity, g);

    // Add edges with capacities
    for (auto& e : edges) {
        auto [u, v, c] = e;
        if (c < 0) throw std::runtime_error("MinCutProblem: negative edges are not allowed!");
        auto e1 = add_edge(u, v, g).first;
        auto e1r = add_edge(v, u, g).first;
        auto e2 = add_edge(v, u, g).first; 
        auto e2r = add_edge(u, v, g).first;
        capacity[e1] = capacity[e2] = c; 
        capacity[e1r] = capacity[e2r] = 0;
        rev[e1] = e1r;
        rev[e1r] = e1;
        rev[e2] = e2r;
        rev[e2r] = e2;
    };

    // Compute max flow
    int minCut = push_relabel_max_flow(g, s, t);

    std::vector<bool> visited(vertices, false);

    std::function<void(int)> dfs = [&](int v) -> void {
        if (visited[v]) return;
        visited[v] = true;
        for (auto [ei, e_end] = out_edges(v, g); ei != e_end; ++ei) {
            if (residual_capacity[*ei] > 0) {
                dfs(target(*ei, g));
            }
        }
    };

    // Run DFS on the residual network to determine the min cut subset
    dfs(s);

    return std::make_pair(minCut, visited);
}


std::pair<int, std::vector<bool>> solveGlobalMinCut(std::vector<std::tuple<int,int,int>> edges) {
    using namespace boost;

    typedef adjacency_list<vecS, vecS, undirectedS, no_property, property<edge_weight_t, int>> Graph;
    typedef property_map<Graph, edge_weight_t>::type WeightMap;

    Graph g;
    for (auto& e : edges) {
        auto [u, v, c] = e;
        if (c < 0) throw std::runtime_error("MinCutProblem: negative edges are not allowed!");
        add_edge(u, v, c, g);
    };

    WeightMap weight_map = get(edge_weight, g);

    std::vector<bool> parity(num_vertices(g));
    auto parity_map = make_iterator_property_map(parity.begin(), get(vertex_index, g));

    int minCut = stoer_wagner_min_cut(g, weight_map, boost::parity_map(parity_map));

    return std::make_pair(minCut, parity);
}




template <typename S = int> // domain S of samples
class UnorderedTriple { // unordered triple (sorted in the ascending order)
    std::array<S, 3> s;

public:

    UnorderedTriple(S s1, S s2, S s3) {
        if (s1 == s2 || s1 == s2 || s2 == s3) throw std::runtime_error("Unordered triple cannot contain the same elements!");
        s = {s1, s2, s3};
        std::sort(std::begin(s), std::end(s));
    }

    S operator[](int index) const { // only getter, not setter (S&)
        return s[index];
    }

    bool operator<(const UnorderedTriple<S>& other) const {
        if (s[0] < other[0]) return true;
        if (s[0] > other[0]) return false;
        if (s[1] < other[1]) return true;
        if (s[1] > other[1]) return false;
        return s[2] < other[2];
    }

    bool contains(S s1) {
        return (s[0] == s1 || s[1] == s1 || s[2] == s1);
    }
};

template <typename S = int> // domain S of samples
class UnorderedPair { // unordered pair (sorted in the ascending order)
    std::array<S, 2> s;

public:

    UnorderedPair(S s1, S s2) {
        if (s1 == s2 || s1 == s2) throw std::runtime_error("Unordered pair cannot contain the same elements!");
        s = {s1, s2};
        std::sort(std::begin(s), std::end(s));
    }

    S operator[](int index) const { // only getter, not setter (S&)
        return s[index];
    }

    bool operator<(const UnorderedPair<S>& other) const {
        if (s[0] < other[0]) return true;
        if (s[0] > other[0]) return false;
        return s[1] < other[1];
    }

    bool contains(S s1) {
        return (s[0] == s1 || s[1] == s1);
    }
};




template <typename S = int> // domain S of samples
class CubicSetPartitionProblem {

    std::vector<S> samples;
    
    int sampleCount; 

    std::vector<std::vector<std::pair<int, int>>> relevantTriples;
    std::map<UnorderedTriple<>, int> tripleCosts;

    std::vector<std::vector<int>> relevantPairs;
    std::map<UnorderedPair<>, int> pairCosts;

    std::vector<int> indexClusterMapping;
    int resultingCost;

    explicit CubicSetPartitionProblem(
        int sampleCount,
        std::vector<std::vector<std::pair<int, int>>> relevantTriples,
        std::map<UnorderedTriple<>, int> tripleCosts,
        std::vector<std::vector<int>> relevantPairs,
        std::map<UnorderedPair<>, int> pairCosts
    ) : relevantTriples(relevantTriples), tripleCosts(tripleCosts),
    relevantPairs(relevantPairs), pairCosts(pairCosts) {
        this->sampleCount = sampleCount;
    }

    std::pair<int, std::vector<bool>> solveMinCutForIndexSubset(bool globalMinCut, std::vector<bool> indexSubset, bool invertCosts, int s = 0, int t = 0) {
        // apply proposition 4.2
        int vertices = std::count(std::begin(indexSubset), std::end(indexSubset), true);
        std::vector<int> indexMapping(indexSubset.size());
        std::vector<int> invIndexMapping(vertices);
        for (int i = 0, sampleNode = 0; i < indexSubset.size(); i++) {
            if (indexSubset[i]) {
                indexMapping[i] = sampleNode;
                invIndexMapping[sampleNode] = i;
                sampleNode++;
            }
        }

        // compute the adjacency matrix by transforming triples (the costs in the matrix are not divided by 2 to avoid floating numbers)
        std::vector<std::vector<int>> adjMatrix(vertices, std::vector(vertices, 0));
        for (int i = 0; i < indexSubset.size(); i++) {
            if (!indexSubset[i]) continue;
            for (int j : relevantPairs[i]) {
                if (i > j) continue; // consider only (i < j)
                if (!indexSubset[j]) continue;
                auto indexPair = UnorderedPair<>(i, j); // sorted indices
                int c = pairCosts[indexPair] * (invertCosts ? -1 : 1);
                int i_node = indexMapping[indexPair[0]]; 
                int j_node = indexMapping[indexPair[1]];
                adjMatrix[i_node][j_node] += c;
            }
            for (auto [j, k] : relevantTriples[i]) {
                if (i > j) continue; // consider only (i < j < k)
                if (!indexSubset[j] || !indexSubset[k]) continue;
                auto indexTriple = UnorderedTriple<>(i, j, k); // sorted indices
                int c = tripleCosts[indexTriple] * (invertCosts ? -1 : 1);
                int i_node = indexMapping[indexTriple[0]]; 
                int j_node = indexMapping[indexTriple[1]];
                int k_node = indexMapping[indexTriple[2]];
                adjMatrix[i_node][j_node] += c;
                adjMatrix[i_node][k_node] += c;
                adjMatrix[j_node][k_node] += c; 
            }
        }

        // create adjacency list from the adjacency matrix
        std::vector<std::tuple<int,int,int>> edges;
        for (int i_node = 0; i_node < vertices; i_node++) {
            for (int j_node = i_node + 1; j_node < vertices; j_node++) {
                int c = adjMatrix[i_node][j_node];
                if (c) {
                    edges.push_back(std::make_tuple(i_node, j_node, c));
                }
            }
        } 

        // solve the MinCut problem
        std::pair<int, std::vector<bool>> solution;
        if (globalMinCut) {
            solution = solveGlobalMinCut(edges);
        } else {
            solution = solveMinCut(vertices, edges, indexMapping[s], indexMapping[t]);
        }
        int minCut = solution.first;
        std::vector<bool> partition = solution.second;

        // transform the results to the original domain
        std::vector<bool> partitionOnSamples(indexSubset.size(), false);
        for (int i_node = 0; i_node < vertices; i_node++) {
            partitionOnSamples[invIndexMapping[i_node]] = partition[i_node];
        }

        // divide the MinCut result by 2 because all triples have been transformed with a double cost
        return std::make_pair(minCut/2, partitionOnSamples);
    }

    CubicSetPartitionProblem<S> createIndependentCutSubproblem(std::vector<int> subsamples) {
        int subSampleCount = subsamples.size();
        
        // subsamples as a set
        std::set<int> subsampleSet;
        for (auto i : subsamples) {
            subsampleSet.insert(i);
        }

        // index of samples in the subproblem
        std::vector<int> indexOf(sampleCount);
        for (int ind = 0; ind < subSampleCount; ind++) {
            indexOf[subsamples[ind]] = ind;
        }

        // filter relevant pairs for th subproblem
        std::vector<std::vector<int>> subRelevantPairs(subSampleCount);
        std::map<UnorderedPair<>, int> subPairCosts;
        for (int i : subsamples) {
            for (int j : relevantPairs[i]) {
                if (subsampleSet.count(j)) {
                    subRelevantPairs[indexOf[i]].push_back(indexOf[j]);
                    subPairCosts[UnorderedPair<>(indexOf[i], indexOf[j])] = pairCosts[UnorderedPair<>(i, j)];
                }
            }
        }
        // filter relevant triples for the subproblem
        std::vector<std::vector<std::pair<int, int>>> subRelevantTriples(subSampleCount);
        std::map<UnorderedTriple<>, int> subTripleCosts;
        for (int i : subsamples) {
            for (auto [j, k] : relevantTriples[i]) {
                if (subsampleSet.count(j) && subsampleSet.count(k)) {
                    subRelevantTriples[indexOf[i]].push_back(std::make_pair(indexOf[j], indexOf[k]));
                    subTripleCosts[UnorderedTriple(indexOf[i], indexOf[j], indexOf[k])] = tripleCosts[UnorderedTriple<>(i, j, k)];
                }
            }
        }
        // create the subproblem
        return CubicSetPartitionProblem(
            subSampleCount,
            subRelevantTriples,
            subTripleCosts,
            subRelevantPairs,
            subPairCosts
        );
    }

    bool applyIndependentSubproblemCut() {
        std::vector<std::vector<int>> partition;
        std::vector<bool> processed(sampleCount, false);
        for (int i = 0; i < sampleCount; i++) {
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
                    UnorderedPair<> indexPair(current, j);
                    int c = pairCosts[indexPair];
                    if (c >= 0) continue; 
                    if (!processed[j]) {
                        processed[j] = true;
                        chosen.push_back(j);
                        q.push(j);
                    }
                }
                for (auto [j, k] : relevantTriples[current]) {
                    UnorderedTriple<> indexTriple(current, j, k);
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

        int clusterOffset = 0;
        for (auto backIndexing : partition) {
            // create and solve the independent subproblems
            auto subproblem = createIndependentCutSubproblem(backIndexing);
            subproblem.solve();
            // accumulate the results
            resultingCost += subproblem.getResultingCost();
            int subClusterCount = 1;
            auto subIndexClusterMapping = subproblem.getIndexClusterMapping();
            for (int ind = 0; ind < subIndexClusterMapping.size(); ind++) {
                int subcluster = subIndexClusterMapping[ind];
                indexClusterMapping[backIndexing[ind]] = subcluster + clusterOffset;
                subClusterCount = std::max(subClusterCount, subcluster + 1);
            }
            clusterOffset += subClusterCount;
        }
        return true;
    }

    void createSolveAccumulateJoinSubproblem(std::vector<int> joinSamples) {
        // aplying proposition 5.1
        int subSampleCount = sampleCount - joinSamples.size() + 1;
        int indexOfJoint = subSampleCount - 1;

        // subsamples as a set
        std::set<int> joinSamplesSet;
        for (auto i : joinSamples) {
            joinSamplesSet.insert(i);
        }

        // index of samples in the subproblem
        std::vector<std::vector<int>> backIndexing(subSampleCount);
        std::vector<int> indexOf(sampleCount);
        for (int i = 0, ind = 0; i < sampleCount; i++) {
            if (joinSamplesSet.count(i)) {
                indexOf[i] = indexOfJoint;
                backIndexing[indexOfJoint].push_back(i);
            } else {
                indexOf[i] = ind;
                backIndexing[ind].push_back(i);
                ind++;
            }
        }

        // filter relevant pairs for the subproblem which have no relations to the samples being joint
        std::vector<std::vector<int>> subRelevantPairs(subSampleCount);
        std::map<UnorderedPair<>, int> subPairCosts;
        for (int i = 0; i < sampleCount; i++) {
            if (joinSamplesSet.count(i)) continue; // skip the relations for the joint set
            for (int j : relevantPairs[i]) {
                if (joinSamplesSet.count(j)) continue; // skip the relations for the joint set
                subRelevantPairs[indexOf[i]].push_back(indexOf[j]);
                subPairCosts[UnorderedPair<>(indexOf[i], indexOf[j])] = pairCosts[UnorderedPair<>(i, j)];
            }
        }
        
        // filter relevant triples for the subproblem which have no relations to the samples being joint
        std::vector<std::vector<std::pair<int, int>>> subRelevantTriples(subSampleCount);
        std::map<UnorderedTriple<>, int> subTripleCosts;
        for (int i = 0; i < sampleCount; i++) {
            if (joinSamplesSet.count(i)) continue; // skip the relations for the joint set
            for (auto [j, k] : relevantTriples[i]) {
                if (joinSamplesSet.count(j) || joinSamplesSet.count(k)) continue; // skip the relations for the joint set
                subRelevantTriples[indexOf[i]].push_back(std::make_pair(indexOf[j], indexOf[k]));
                subTripleCosts[UnorderedTriple(indexOf[i], indexOf[j], indexOf[k])] = tripleCosts[UnorderedTriple<>(i, j, k)];
            }
        }

        // compute the costs to the joint subset as well as the inner joining cost
        for (int i : joinSamples) {
            for (int j : relevantPairs[i]) {
                if (joinSamplesSet.count(j)) {
                    if (i < j) resultingCost += pairCosts[UnorderedPair<>(i, j)]; // consider one direction (i, j) and skip (j, i)
                } else {
                    UnorderedPair<> indexPair(indexOf[i], indexOf[j]);
                    if (!subPairCosts.count(indexPair)) subPairCosts[indexPair] = 0;
                    subRelevantPairs[indexOf[i]].push_back(indexOf[j]);
                    subRelevantPairs[indexOf[j]].push_back(indexOf[i]);
                    subPairCosts[indexPair] += pairCosts[UnorderedPair<>(i, j)];
                }
            }
        }
        for (int i : joinSamples) {
            for (auto [j, k] : relevantTriples[i]) {
                bool innerJ = joinSamplesSet.count(j);
                bool innerK = joinSamplesSet.count(k);
                if (innerJ && innerK) {
                    // (i < j < k) consider only (i, j, k) and skip the others
                    if (i < j) resultingCost += tripleCosts[UnorderedTriple<>(i, j, k)];
                } else if (innerJ || innerK) {
                    int outer, inner;
                    if (innerJ) {
                        inner = j;
                        outer = k;
                    }
                    if (innerK) {
                        inner = k;
                        outer = j;
                    }
                    if (i > inner) continue; // consider one direction (i, inner) and skip (inner, i)
                    UnorderedPair<> indexPair(indexOf[i], indexOf[outer]);
                    if (!subPairCosts.count(indexPair)) subPairCosts[indexPair] = 0;
                    subRelevantPairs[indexOf[i]].push_back(indexOf[outer]);
                    subRelevantPairs[indexOf[outer]].push_back(indexOf[i]);
                    subPairCosts[indexPair] += tripleCosts[UnorderedTriple<>(i, j, k)];
                } else { // j and k are outer
                    UnorderedTriple<> indexTriple(indexOf[i], indexOf[j], indexOf[k]);
                    if (!subTripleCosts.count(indexTriple)) subTripleCosts[indexTriple] = 0;
                    subRelevantTriples[indexOf[i]].push_back(std::make_pair(indexOf[j], indexOf[k]));
                    subRelevantTriples[indexOf[j]].push_back(std::make_pair(indexOf[k], indexOf[i])); // indexOf[i] > indexOf[k] by definition above
                    subRelevantTriples[indexOf[k]].push_back(std::make_pair(indexOf[j], indexOf[i])); // indexOf[i] > indexOff[j] by definition above
                    subTripleCosts[indexTriple] += tripleCosts[UnorderedTriple<>(i, j, k)];
                }
            } 
        }

        // create and solve the subproblem
        auto subproblem = CubicSetPartitionProblem(
            subSampleCount,
            subRelevantTriples,
            subTripleCosts,
            subRelevantPairs,
            subPairCosts
        );
        subproblem.solve();

        // accumulate the results (resulting cost for the join has been accumulated)
        resultingCost += subproblem.getResultingCost();
        auto subIndexClusterMapping = subproblem.getIndexClusterMapping();
        for (int ind = 0; ind < subSampleCount; ind++) {
            int subcluster = subIndexClusterMapping[ind];
            for (int originalIndex : backIndexing[ind]) {
                indexClusterMapping[originalIndex] = subcluster;
            }
        }
    }

    bool applyBipartiteSubsetJoinForIndexSubset(std::vector<bool> &indexSubset) {   
        // compute rhs
        int lhsLowerBound = 0; // avoid MinCut computation for lhs>rhs and the edge cases where lhs=0 (rhs<0 because of 3.1) or avoid
        int singleR = 0, doubleR = 0;
        for (int i = 0; i < sampleCount; i++) {
            if (!indexSubset[i]) continue; // i is in R
            for (auto j : relevantPairs[i]) {
                int c = pairCosts[UnorderedPair<>(i, j)];
                if (indexSubset[j]) {
                    if (i < j) lhsLowerBound += c;
                } else if (c < 0) {
                    singleR += c;
                }
            }
            for (auto [j, k] : relevantTriples[i]) {
                int c = tripleCosts[UnorderedTriple<>(i, j, k)];
                if (indexSubset[j] && indexSubset[k]) {
                    if (i < j) lhsLowerBound += c;
                    continue; // j or k must be not in R
                }
                if (c > 0) continue; // omit positive costs
                if (!indexSubset[k] && !indexSubset[j]) {
                    singleR += c;
                } else {
                    doubleR += c; // since two elements of the triple are in R, the cost will be added twice
                }
            }
        }
        int rhs = singleR + doubleR/2;
        if (lhsLowerBound > rhs) return false;

        auto [minCut, partition] = solveMinCutForIndexSubset(true, indexSubset, true);
        int lhs = -minCut; 
        if (lhs <= rhs) {
            std::vector<int> joinSamples;
            for (int i = 0; i < sampleCount; i++) {
                if (indexSubset[i]) joinSamples.push_back(i);
            }
            createSolveAccumulateJoinSubproblem(joinSamples);
            return true;
        }  
        return false;
    }

    bool applyBipartiteSubsetJoin() {
        // subset join criterion 3.11
        // heuristically construct and check the candidate sets R for possible joining
        for (int i = 0; i < sampleCount; i++) {
            for (int j = i + 1; j < sampleCount; j++) {
                if (pairCosts.count(UnorderedPair<>(i, j)) && pairCosts[UnorderedPair<>(i, j)] > 0) continue;

                std::vector<bool> subsetR(sampleCount, false); // R doesn't have to be a connected component
                subsetR[i] = subsetR[j] = true;
                if (applyBipartiteSubsetJoinForIndexSubset(subsetR)) return true;

                std::set<int> candidates;
                for (int k = 0; k < sampleCount; k++) {
                    if (!subsetR[k]) {
                        candidates.insert(k);
                    }
                }

                std::function<int(int)> computeOffset = [&](int i) {
                    // positive offset if not mergeable
                    if (subsetR[i]) return 1;
                    int offset = 0;
                    for (auto j : relevantPairs[i]) {
                        if (!subsetR[j]) continue;
                        int c = pairCosts[UnorderedPair<>(i, j)];
                        if (c > 0) {
                            return 1;
                        } else {
                            offset += c;
                        }
                    }
                    for (auto [j, k] : relevantTriples[i]) {
                        if (!subsetR[j] || !subsetR[k]) continue; // 2 of 3 triple elements are already in R
                        int c = tripleCosts[UnorderedTriple<>(i, j, k)];
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
                        subsetR[bestCandidate] = true;
                        candidates.erase(bestCandidate);
                        if (applyBipartiteSubsetJoinForIndexSubset(subsetR)) return true;
                    }
                }
            }
        }
        return false;
    }
    
    
public: 
    explicit CubicSetPartitionProblem(const std::vector<S>& givenSamples, const std::function<int(UnorderedTriple<S>)> costCB)
    : samples(givenSamples) {
        sampleCount = samples.size();
        // init relevant pairs
        relevantPairs.resize(sampleCount, {});
        // init relevant triples
        relevantTriples.resize(sampleCount, {});
        for (int i = 0; i < sampleCount; i++) {
            for (int j = i + 1; j < sampleCount; j++) {
                for (int k = j + 1; k < sampleCount; k++) {
                    UnorderedTriple<> indexTriple(i, j, k);
                    UnorderedTriple<S> sampleTriple(samples[i], samples[j], samples[k]);
                    int c = costCB(sampleTriple);
                    if (c) {
                        tripleCosts[indexTriple] = c;
                        relevantTriples[i].push_back(std::make_pair(j, k));
                        relevantTriples[j].push_back(std::make_pair(i, k));
                        relevantTriples[k].push_back(std::make_pair(i, j));
                    }
                }
            }
        }
    }

    std::map<S, int> getClusterMapping() {
        if (samples.empty()) throw std::runtime_error("Cluster mapping for samples is available only for the original problem, not subproblems!");
        std::map<S, int> result;
        for (int i = 0; i < sampleCount; i++) {
            result[samples[i]] = indexClusterMapping[i];
        }
        return result;
    } 

    std::vector<int> getIndexClusterMapping() {
        return indexClusterMapping;
    } 

    int getResultingCost() {
        return resultingCost;
    }

    void solve() {
        if (!sampleCount) throw std::runtime_error("Cannot solve a cubic set partition problem with no samples!");

        // init the problem results (all elements are joint by default)
        indexClusterMapping.resize(sampleCount, 0);
        resultingCost = 0;
        
        if (sampleCount == 1) return; // trivial problem

        // apply partial optimality conditions (solve the subproblems)
        if (applyIndependentSubproblemCut()) return;
        if (applyBipartiteSubsetJoin()) return;
        // TODO: apply other conditions if-return
        
        // current problem could not be reduced to subproblems
        std::cout << "WARNING: found a non-trivial problem (2+ samples) that has been reduced to subproblems!" << std::endl;
        return;
    }

};



// int cost(UnorderedTriple<char> t) {
//     if (t[0] == 'a' && t[1] == 'b' && t[2] == 'c')
//         return -5;
//     if (t[0] == 'b' && t[1] == 'c' && t[2] == 'd')
//         return -7;
//     if (t[0] == 'a' && t[1] == 'b' && t[2] == 'd')
//         return 10;
//     return 0;
// }

int cost(UnorderedTriple<char> t) {
    if (t[0] == 'a' && t[1] == 'b' && t[2] == 'c') return -1;
    if (t[0] == 'a' && t[1] == 'c' && t[2] == 'd') return -15;
    if (t[0] == 'd' && t[1] == 'e' && t[2] == 'h') return 50;
    if (t[0] == 'e' && t[1] == 'f' && t[2] == 'h') return -50;
    if (t[0] == 'f' && t[1] == 'g' && t[2] == 'i') return -30;
    if (t[0] == 'd' && t[1] == 'h' && t[2] == 'k') return -2;
    if (t[0] == 'i' && t[1] == 'k' && t[2] == 'l') return -4;
    if (t[0] == 'j' && t[1] == 'l' && t[2] == 'm') return -10;
    return 0;
}


int main() {
    std::vector<char> samples = {'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm'};
    CubicSetPartitionProblem<char> problem(samples, cost);
    problem.solve();
    for (auto [sample, cluster] : problem.getClusterMapping()) {
        std::cout << sample << " -> " << cluster << std::endl;
    }
    std::cout << problem.getResultingCost() << std::endl;

    return 0;
}