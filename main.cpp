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


int solveMinCut(int vertices, std::vector<std::tuple<int,int,int>> edges, int s, int t) {
    // handled as bigraph!
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

    return minCut;
}


int solveGlobalMinCut(std::vector<std::tuple<int,int,int>> edges) {
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

    return minCut;
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

    std::vector<std::vector<bool>> labelFixed, labelValue;
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

    int solveMinCutForIndexSubset(std::vector<bool> indexSubset, bool takeNegativeCosts, bool takePositiveCosts,  bool globalMinCut, int source = 0, std::vector<int> sinks = std::vector<int>({0})) {
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
                auto indexPair = UnorderedPair<>(i, j); // sorted indices
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
                auto indexTriple = UnorderedTriple<>(i, j, k); // sorted indices
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
            minCut = solveGlobalMinCut(edges);
        } else {
            int superSink = vertices;
            for (auto sink : sinks) {
                edges.push_back(std::make_tuple(indexMapping[sink], superSink, sumCosts));
            }
            minCut = solveMinCut(vertices + 1, edges, indexMapping[source], superSink);
        }

        // divide the MinCut result by 2 because all triples have been transformed with a double cost
        return minCut/2;
    }

    CubicSetPartitionProblem<S> createIndependentCutSubproblem(std::vector<int> subsamples) {
        int subSampleCount = subsamples.size();
        
        // subsamples as a set
        std::vector<bool> subsetR(sampleCount, false);
        for (auto i : subsamples) {
            subsetR[i] = true;
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
                if (subsetR[j]) {
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
                if (subsetR[j] && subsetR[k]) {
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
        std::cout << "Applying the independent subproblem cut (proposition 3.1)" << std::endl;

        // fix the labels
        for (int subsetI = 0; subsetI < partition.size(); subsetI++) {
            for (int subsetJ = subsetI + 1; subsetJ < partition.size(); subsetJ++) {
                for (auto i : partition[subsetI]) {
                    for (auto j : partition[subsetJ]) {
                        labelFixed[i][j] = labelFixed[j][i] = true;
                        labelValue[i][j] = labelValue[j][i] = false;
                    }
                }
            }
        }

        // create, solve the subproblems, accumulate the results
        int clusterOffset = 0;
        for (auto backIndexing : partition) {
            // create and solve the independent subproblems
            auto subproblem = createIndependentCutSubproblem(backIndexing);
            subproblem.solve();
            // accumulate the results
            resultingCost += subproblem.getResultingCost();
            // accumulate the labels
            auto subLabelFixed = subproblem.getIndexLabelFixed();
            auto subLabelValue = subproblem.getIndexLabelValue();
            for (int i = 0; i < subLabelFixed.size(); i++) {
                for (int j = i + 1; j < subLabelFixed.size(); j++) {
                    int originalI = backIndexing[i], originalJ = backIndexing[j];
                    labelFixed[originalI][originalJ] = labelFixed[originalJ][originalI] = subLabelFixed[i][j];
                    labelValue[originalI][originalJ] = labelValue[originalJ][originalI] = subLabelValue[i][j];
                }
            } 
        }
        
        return true;
    }

    void createSolveAccumulateJoinSubproblem(std::vector<bool> subsetR) {
        std::vector<int> joinSamples;
        for (int i = 0; i < sampleCount; i++) {
            if (subsetR[i]) joinSamples.push_back(i);
        }

        // aplying proposition 5.1
        int subSampleCount = sampleCount - joinSamples.size() + 1;
        int indexOfJoint = subSampleCount - 1;

        // index of samples in the subproblem
        std::vector<std::vector<int>> backIndexing(subSampleCount);
        std::vector<int> indexOf(sampleCount);
        for (int i = 0, ind = 0; i < sampleCount; i++) {
            if (subsetR[i]) {
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
            if (subsetR[i]) continue; // skip the relations for the joint set
            for (int j : relevantPairs[i]) {
                if (subsetR[j]) continue; // skip the relations for the joint set
                subRelevantPairs[indexOf[i]].push_back(indexOf[j]);
                subPairCosts[UnorderedPair<>(indexOf[i], indexOf[j])] = pairCosts[UnorderedPair<>(i, j)];
            }
        }
        
        // filter relevant triples for the subproblem which have no relations to the samples being joint
        std::vector<std::vector<std::pair<int, int>>> subRelevantTriples(subSampleCount);
        std::map<UnorderedTriple<>, int> subTripleCosts;
        for (int i = 0; i < sampleCount; i++) {
            if (subsetR[i]) continue; // skip the relations for the joint set
            for (auto [j, k] : relevantTriples[i]) {
                if (subsetR[j] || subsetR[k]) continue; // skip the relations for the joint set
                subRelevantTriples[indexOf[i]].push_back(std::make_pair(indexOf[j], indexOf[k]));
                subTripleCosts[UnorderedTriple(indexOf[i], indexOf[j], indexOf[k])] = tripleCosts[UnorderedTriple<>(i, j, k)];
            }
        }

        // compute the costs to the joint subset as well as the inner joining cost
        for (int i : joinSamples) {
            for (int j : relevantPairs[i]) {
                if (subsetR[j]) {
                    if (i < j) resultingCost += pairCosts[UnorderedPair<>(i, j)]; // consider one direction (i, j) and skip (j, i)
                } else {
                    UnorderedPair<> indexPair(indexOf[i], indexOf[j]);
                    if (!subPairCosts.count(indexPair)) {
                        subPairCosts[indexPair] = 0;
                        subRelevantPairs[indexOf[i]].push_back(indexOf[j]);
                        subRelevantPairs[indexOf[j]].push_back(indexOf[i]);
                    }
                    subPairCosts[indexPair] += pairCosts[UnorderedPair<>(i, j)];
                }
            }
        }
        for (int i : joinSamples) {
            for (auto [j, k] : relevantTriples[i]) {
                bool innerJ = subsetR[j], innerK = subsetR[k];
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
                    if (!subPairCosts.count(indexPair)) {
                        subPairCosts[indexPair] = 0;
                        subRelevantPairs[indexOf[i]].push_back(indexOf[outer]);
                        subRelevantPairs[indexOf[outer]].push_back(indexOf[i]);
                    }
                    subPairCosts[indexPair] += tripleCosts[UnorderedTriple<>(i, j, k)];
                } else { // j and k are outer
                    UnorderedTriple<> indexTriple(indexOf[i], indexOf[j], indexOf[k]);
                    if (!subTripleCosts.count(indexTriple)) {
                        subTripleCosts[indexTriple] = 0;
                        subRelevantTriples[indexOf[i]].push_back(std::make_pair(indexOf[j], indexOf[k]));
                        subRelevantTriples[indexOf[j]].push_back(std::make_pair(indexOf[k], indexOf[i])); // indexOf[i] > indexOf[k] by definition above
                        subRelevantTriples[indexOf[k]].push_back(std::make_pair(indexOf[j], indexOf[i])); // indexOf[i] > indexOff[j] by definition above
                    }
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
        // fix the labels for join
        for (int indI = 0; indI < joinSamples.size(); indI++) {
            for (int indJ = indI + 1; indJ < joinSamples.size(); indJ++) {
                int i = joinSamples[indI], j = joinSamples[indJ];
                labelFixed[i][j] = labelFixed[j][i] = true;
                labelValue[i][j] = labelValue[j][i] = true;
            }
        }
        // accumulate the labels
        auto subLabelFixed = subproblem.getIndexLabelFixed();
        auto subLabelValue = subproblem.getIndexLabelValue();
        for (int i = 0; i < subLabelFixed.size(); i++) {
            for (int j = i + 1; j < subLabelFixed.size(); j++) {
                for (int originalI : backIndexing[i]) {
                    for (int originalJ : backIndexing[j]) {
                        labelFixed[originalI][originalJ] = labelFixed[originalJ][originalI] = subLabelFixed[i][j];
                        labelValue[originalI][originalJ] = labelValue[originalJ][originalI] = subLabelValue[i][j];
                    }
                }   
            }
        } 
    }

    bool checkSubsetJoinForIndexSubset(std::vector<bool> &indexSubset) {   
        // compute rhs
        int lhsLowerBound = 0; // avoid MinCut computation for lhs>rhs (in particular the edge cases with lhs=0, rhs<0 because of 3.1 applied before)
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

        int lhs = -solveMinCutForIndexSubset(indexSubset, true, false, true);
        return (lhs <= rhs);
    }

    bool applySubsetJoin() {
        // subset join criterion 3.11
        // heuristically construct and check the candidate sets R for possible joining
        for (int i = 0; i < sampleCount; i++) {
            for (int j = i + 1; j < sampleCount; j++) {
                if (pairCosts.count(UnorderedPair<>(i, j)) && pairCosts[UnorderedPair<>(i, j)] > 0) continue;

                std::vector<bool> subsetR(sampleCount, false); // R doesn't have to be a connected component
                subsetR[i] = subsetR[j] = true;
                std::vector<bool> joinR;
                if (checkSubsetJoinForIndexSubset(subsetR)) joinR = subsetR;

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
                        if (checkSubsetJoinForIndexSubset(subsetR)) joinR = subsetR;
                    }
                }
                if (!joinR.empty()) {
                    std::cout << "Applying the subset join (proposition 3.11)" << std::endl;
                    createSolveAccumulateJoinSubproblem(joinR);
                    return true;
                }
            }
        }
        return false;
    }
    
    bool applyPairJoin() {
        // proposition 3.4
        for (int i = 0; i < sampleCount; i++) {
            for (int j = i + 1; j < sampleCount; j++) {
                // compute lhs
                int lhs = 0;
                UnorderedPair<> indexPair(i, j);
                if (pairCosts.count(indexPair)) {
                    int c = pairCosts[indexPair];
                    if (c < 0) lhs += -2*c;
                }
                for (auto [k1, k2] : relevantTriples[i]) {
                    if (k1 != j && k2 != j) continue; // k1 or k2 is j
                    int c = tripleCosts[UnorderedTriple<> (i, k1, k2)];
                    if (c < 0) lhs += -c;
                }
                // compute rhs (assume that the underlying graph is connected after applying the proposition 3.1)
                int rhs = solveMinCutForIndexSubset(std::vector<bool>(sampleCount, true), true, true, false, i, {j});

                if (lhs >= rhs) {
                    std::vector<bool> subsetR(sampleCount, false);
                    subsetR[i] = subsetR[j] = true;
                    std::cout << "Applying the pair join (proposition 3.4)" << std::endl;
                    // std::cout << i << " " << j << std::endl;
                    // std::cout << lhs << " vs " << rhs << std::endl;
                    createSolveAccumulateJoinSubproblem(subsetR);
                    return true;
                }
            }
        }
        return false;
    }

    bool applyComplexPairJoin() {
        // proposition 3.6
        for (int i = 0; i < sampleCount; i++) {
            for (int k = i + 1; k < sampleCount; k++) {
                for (int j = 0; j < sampleCount; j++) {
                    if (j == i || j == k) continue;
                    // compute lhs1, lhs2, lhs3
                    int lhs1 = 0, lhs2 = 0, lhs3 = 0;
                    UnorderedTriple<> indexTriple(i, j, k);
                    if (tripleCosts.count(indexTriple)) {
                        int c = tripleCosts[indexTriple];
                        lhs3 += c;
                        if (c < 0) {
                            lhs1 += -c;
                            lhs2 += -c; 
                        }
                    }
                    UnorderedPair<> indexPairIJ(i, j), indexPairIK(i, k), indexPairJK(j, k);
                    if (pairCosts.count(indexPairIJ)) {
                        int c = pairCosts[indexPairIJ];
                        lhs3 += c;
                        if (c < 0) lhs1 += -2*c; 
                    }
                    if (pairCosts.count(indexPairJK)) {
                        int c = pairCosts[indexPairJK];
                        lhs3 += c;
                        if (c < 0) lhs2 += -2*c; 
                    }
                    if (pairCosts.count(indexPairIK)) {
                        int c = pairCosts[indexPairIK];
                        lhs3 += c;
                        if (c < 0) {
                            lhs1 += -2*c;
                            lhs2 += -2*c;
                        }
                    }
                    for (auto [p, q] : relevantTriples[i]) {
                        if (p == j || p == k || q == j || q == k) {
                            int c = tripleCosts[UnorderedTriple<>(i, p, q)];
                            if (c < 0) lhs1 += -c;
                        }
                    }
                    for (auto [p, q] : relevantTriples[k]) {
                        if (p == j || p == i || q == j || q == i) {
                            int c = tripleCosts[UnorderedTriple<>(k, p, q)];
                            if (c < 0) lhs2 += -c;
                        }
                    }
                    // compute rhs3
                    int rhs3 = 0;
                    for (auto [p, q] : relevantTriples[i]) {
                        if (p == j || p == k || q == j || q == k) continue;
                        int c = tripleCosts[UnorderedTriple<>(i, p, q)];
                        rhs3 -= abs(c);
                    }
                    for (auto [p, q] : relevantTriples[j]) {
                        if (p == i || p == k || q == i || q == k) continue;
                        int c = tripleCosts[UnorderedTriple<>(j, p, q)];
                        rhs3 -= abs(c);
                    }
                    for (auto [p, q] : relevantTriples[k]) {
                        if (p == i || p == j || q == i || q == j) continue;
                        int c = tripleCosts[UnorderedTriple<>(k, p, q)];
                        rhs3 -= abs(c);
                    }
                    for (auto p : relevantPairs[i]) {
                        if (p == j || p == k) continue;
                        int c = pairCosts[UnorderedPair(i, p)];
                        rhs3 -= abs(c);
                    }
                    for (auto p : relevantPairs[j]) {
                        if (p == i || p == k) continue;
                        int c = pairCosts[UnorderedPair(j, p)];
                        rhs3 -= abs(c);
                    }
                    for (auto p : relevantPairs[k]) {
                        if (p == i || p == j) continue;
                        int c = pairCosts[UnorderedPair(k, p)];
                        rhs3 -= abs(c);
                    }
                    // check the 3d condition
                    if (!(lhs3 <= rhs3)) continue; 
                    // check the 1st and the 2d condition
                    int rhs1 = solveMinCutForIndexSubset(std::vector<bool>(sampleCount, true), true, true, false, i, {j, k});
                    int rhs2 = solveMinCutForIndexSubset(std::vector<bool>(sampleCount, true), true, true, false, k, {i, j});
                    if (lhs1 >= rhs1 && lhs2 >= rhs2) {
                        std::vector<bool> subsetR(sampleCount, false);
                        subsetR[i] = subsetR[k] = true; // join i and k
                        std::cout << "Applying the complex pair join (proposition 3.6)" << std::endl;
                        // std::cout << i << " " << j << " " << k << std::endl;
                        // std::cout << "(" << lhs1 << "," << lhs2 << "," << lhs3 << ") vs ";
                        // std:: cout << "(" << rhs1 << "," << rhs2 << "," << rhs3 << ")" << std::endl;
                        createSolveAccumulateJoinSubproblem(subsetR);
                        return true;
                    }
                }
            }
        }
        return false;
    }

    bool applyExplicitPairJoin() {
        // proposition 3.8
        for (int i = 0; i < sampleCount; i++) {
            for (int j = i + 1; j < sampleCount; j++) {
                int lhs = 0;
                UnorderedPair indexPair(i, j);
                if (pairCosts.count(indexPair)) lhs += pairCosts[indexPair];
                // compute rhs
                int rhs = 0;
                for (int k : relevantPairs[i]) {
                    if (k == j) continue;
                    int c = pairCosts[UnorderedPair<>(i, k)];
                    if (c < 0) rhs += c;
                }
                for (int k : relevantPairs[j]) {
                    if (k == i) continue;
                    int c = pairCosts[UnorderedPair<>(j, k)];
                    if (c < 0) rhs += c;
                }
                for (auto [k1, k2] : relevantTriples[i]) {
                    int c = tripleCosts[UnorderedTriple<>(i, k1, k2)];
                    if (c < 0) rhs += c;
                }
                for (auto [k1, k2] : relevantTriples[j]) {
                    if (k1 == i || k2 == i) continue; // the triples (i, j, *) have already been considered in the previos for-loop
                    int c = tripleCosts[UnorderedTriple<>(j, k1, k2)];
                    if (c < 0) rhs += c;
                }
                // check the condition
                if (lhs <= rhs) {
                    std::vector<bool> subsetR(sampleCount, false);
                    subsetR[i] = subsetR[j] = true;
                    std::cout << "Applying the explicit pair join (proposition 3.8)" << std::endl;
                    std::cout << i << " " << j << std::endl;
                    std::cout << lhs << " vs " << rhs << std::endl;
                    createSolveAccumulateJoinSubproblem(subsetR);
                    return true;
                }
            }
        }
        return false;
    }

    bool applyExplicitPairJoinViaTriple() {
        // proposition 3.9
        for (int i = 0; i < sampleCount; i++) {
            for (int k = i + 1; k < sampleCount; k++) {
                for (int j = 0; j < sampleCount; j++) {
                    if (j == i || j == k) continue;
                    // compute costs for the triple
                    UnorderedPair<> indexPairIJ(i, j), indexPairIK(i, k), indexPairJK(j, k);
                    UnorderedTriple<> indexTriple(i, j, k);
                    int cIJ = 0, cIK = 0, cJK = 0, cIJK = 0;
                    if (pairCosts.count(indexPairIJ)) cIJ = pairCosts[indexPairIJ];
                    if (pairCosts.count(indexPairIK)) cIK = pairCosts[indexPairIK];
                    if (pairCosts.count(indexPairJK)) cJK = pairCosts[indexPairJK];
                    if (tripleCosts.count(indexTriple)) cIJK = tripleCosts[indexTriple];
                    // compute rhs
                    int singleR = 0, doubleR = 0;
                    std::vector<int> elementsR = {i, j, k};
                    std::vector<bool> subsetR(sampleCount, false);
                    subsetR[i] = subsetR[j] = subsetR[k] = true;
                    for (int i1 : elementsR) {
                        for (int j1 : relevantPairs[i1]) {
                            if (subsetR[j1]) continue;
                            int c = pairCosts[UnorderedPair<>(i1, j1)];
                            if (c < 0) singleR += c;
                        }
                        for (auto [j1, k1] : relevantTriples[i1]) {
                            if (subsetR[j1] && subsetR[k1]) continue;
                            int c = tripleCosts[UnorderedTriple<>(i1, j1, k1)];
                            if (c > 0) continue;
                            if (subsetR[j1] xor subsetR[k1]) {
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
                        std::vector<bool> subsetR(sampleCount, false);
                        subsetR[i] = subsetR[k] = true; // join i and k
                        std::cout << "Applying the complex pair join (proposition 3.9)" << std::endl;
                        std::cout << i << " " << j << " " << k << std::endl;
                        std::cout << rhs << std::endl;
                        createSolveAccumulateJoinSubproblem(subsetR);
                        return true;
                    }
                }
            }
        }
        return false;
    }

    bool applyTripleJoin() {
        // proposition 3.5
        // compute lhs base value
        int lhsBase = 0;
        for (int i = 0; i < sampleCount; i++) {
            for (int j : relevantPairs[i]) {
                if (i > j) continue;
                int c = pairCosts[UnorderedPair<>(i, j)];
                if (c > 0) lhsBase -= c;
            }
            for (auto [j, k] : relevantTriples[i]) {
                if (i > j) continue;
                int c = tripleCosts[UnorderedTriple<>(i, j, k)];
                if (c > 0) lhsBase -= c;
            }
        }
        // iterate over all unordered triples (i, j, k)
        for (int i = 0; i < sampleCount; i++) {
            for (int j = i + 1; j < sampleCount; j++) {
                for (int k = j + 1; k < sampleCount; k++) {
                    // compute lhs
                    int lhs = lhsBase;
                    UnorderedPair<> indexPairIJ(i, j), indexPairIK(i, k), indexPairJK(j, k);
                    UnorderedTriple<> indexTriple(i, j, k);
                    int cIJ = 0, cIK = 0, cJK = 0, cIJK = 0;
                    if (pairCosts.count(indexPairIJ)) cIJ = pairCosts[indexPairIJ];
                    if (cIJ < 0) lhs += -2*cIJ;
                    if (pairCosts.count(indexPairIK)) cIK = pairCosts[indexPairIK];
                    if (cIK < 0) lhs += -2*cIK;
                    if (pairCosts.count(indexPairJK)) cJK = pairCosts[indexPairJK];
                    if (cJK < 0) lhs += -cJK;
                    if (tripleCosts.count(indexTriple)) cIJK = tripleCosts[indexTriple];
                    if (cIJK < 0) lhs += -2*cIJK;
                    lhs += std::min(std::min(0, cIJ), std::min(cIK, cJK)); // min for 4 cases
                    // compute rhs
                    int rhs = solveMinCutForIndexSubset(std::vector<bool>(sampleCount, true), true, false, false, i, {j, k});
                    if (lhs >= rhs) {
                        std::vector<bool> subsetR(sampleCount, false);
                        subsetR[i] = subsetR[j] = subsetR[k] = true; // join ijk
                        std::cout << "Applying the triple join (proposition 3.5)" << std::endl;
                        std::cout << i << " " << j << " " << k << std::endl;
                        std::cout << lhs << " vs " << rhs << std::endl;
                        createSolveAccumulateJoinSubproblem(subsetR);
                        return true;
                    }
                }
            }
        }
        return false;
    }

    void applyPairCuts() {
        // proposition 3.2
        for (int i = 0; i < sampleCount; i++) {
            for (int j = i + 1; j < sampleCount; j++) {
                // check if already cut
                if (labelFixed[i][j]) {
                    if (labelValue[i][j]) throw std::runtime_error("Error in the program detected: cannot cut the joint elements! ");
                    else continue;
                }
                // check the condition
                int lhs = 0;
                UnorderedPair<> indexPair(i, j);
                if (pairCosts.count(indexPair)) {
                    int c = pairCosts[indexPair];
                    if (c > 0) lhs += c;
                }
                int rhs = solveMinCutForIndexSubset(std::vector<bool>(sampleCount, true), true, false, false, i, {j});
                if (lhs >= rhs) {
                    std::cout << "Applying the pair cut (proposition 3.2)" << std::endl;
                    std::cout << i << " " << j << std::endl;
                    std::cout << lhs << " vs " << rhs << std::endl;
                    // assume no joint labels in this problems, because it would be solved by a join-subproblem!
                    labelFixed[i][j] = labelFixed[j][i] = true;
                    labelValue[i][j] = labelValue[j][i] = false;
                }
            }
        }
    }

    void applyTripleCuts() {
        // TODO
    }


public: 
    explicit CubicSetPartitionProblem(const std::vector<S>& givenSamples, const std::function<int(UnorderedTriple<S>)> tripleCostCB, const std::function<int(UnorderedPair<S>)> pairCostCB  = [](UnorderedPair<S> p)->int{return 0;})
    : samples(givenSamples) {
        sampleCount = samples.size();
        // init relevant pairs
        relevantPairs.resize(sampleCount, {});
        for (int i = 0; i < sampleCount; i++) {
            for (int j = i + 1; j < sampleCount; j++) {
                UnorderedPair<> indexPair(i, j);
                UnorderedPair<S> samplePair(samples[i], samples[j]);
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
                    UnorderedTriple<> indexTriple(i, j, k);
                    UnorderedTriple<S> sampleTriple(samples[i], samples[j], samples[k]);
                    int c = tripleCostCB(sampleTriple);
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

    std::vector<std::vector<bool>> getIndexLabelValue() {
        return labelValue;
    }

    std::vector<std::vector<bool>> getIndexLabelFixed() {
        return labelFixed;
    }

    bool isSolvedCompletely() {
        for (int i = 0; i < sampleCount; i++) {
            for (int j = i + 1; j < sampleCount; j++) {
                if (!labelFixed[i][j]) return false;
            }
        }
        return true;
    }

    int getResultingCost() {
        return resultingCost;
    }

    void printResultingLabels() {
        std::cout << "  ";
        for (auto s : samples) {
            std::cout << s << " ";
        }
        std::cout << std::endl;
        for (int i = 0; i < samples.size(); i++) {
            std::cout << samples[i] << " ";
            for (int j = 0; j < samples.size(); j++) {
                if (i == j) {
                    std::cout << "-";
                } else if (labelFixed[i][j]) {
                    std::cout << int(labelValue[i][j]);
                } else {
                    std::cout << "x";
                }
                std::cout << " ";
            }
            std::cout << std::endl;
        }
    }
    
    void printResults() {
        bool completeSolution = isSolvedCompletely();
        std::cout << std::endl;
        std::cout << "Problem solved: " << ((completeSolution) ? "completely" : "partially") << std::endl; 
        std::cout << "Resulting cost: " << resultingCost << ((completeSolution) ? "" : " (ignoring the unsolved subproblems)") << std::endl;
        printResultingLabels();
    }

    void solve() {
        if (!sampleCount) throw std::runtime_error("Cannot solve a cubic set partition problem with no samples!");

        // init the problem results
        labelFixed.resize(sampleCount, std::vector<bool>(sampleCount, false));
        labelValue.resize(sampleCount, std::vector<bool>(sampleCount, false));
        resultingCost = 0;

        
        if (sampleCount == 1) return; // trivial problem

        // apply partial optimality conditions (solve the subproblems if needed)
        // if (applyIndependentSubproblemCut()) return;
        // if (applySubsetJoin()) return;
        // if (applyPairJoin()) return;
        // if (applyComplexPairJoin()) return;
        // if (applyExplicitPairJoin()) return;
        // if (applyExplicitPairJoinViaTriple()) return;
        if (applyTripleJoin()) return;
        applyPairCuts();
        applyTripleCuts();
        
        if (!isSolvedCompletely()) {
            // current problem could not be reduced to subproblems
            std::cout << "WARNING: found a non-trivial problem (" << sampleCount << " samples) that has not been reduced to subproblems! " << std::endl;
        }
        return;
    }
};



int cost(UnorderedTriple<char> t) {
    // if (t[0] == 'a' && t[1] == 'b' && t[2] == 'c')
    //     return -5;
    if (t[0] == 'b' && t[1] == 'c' && t[2] == 'd')
        return -22; // -7
    if (t[0] == 'a' && t[1] == 'b' && t[2] == 'd')
        return 10;
    if (t[0] == 'a' && t[1] == 'b' && t[2] == 'e')
        return 10;
    return 0;
}

// int cost(UnorderedTriple<char> t) {
//     // example: 3.1 and 3.11 are sufficient
//     if (t[0] == 'a' && t[1] == 'b' && t[2] == 'c') return -1;
//     if (t[0] == 'a' && t[1] == 'c' && t[2] == 'd') return -15;
//     if (t[0] == 'd' && t[1] == 'e' && t[2] == 'h') return 50;
//     if (t[0] == 'e' && t[1] == 'f' && t[2] == 'h') return -50;
//     if (t[0] == 'f' && t[1] == 'g' && t[2] == 'i') return -30;
//     if (t[0] == 'd' && t[1] == 'h' && t[2] == 'k') return -2;
//     if (t[0] == 'i' && t[1] == 'k' && t[2] == 'l') return -4;
//     if (t[0] == 'j' && t[1] == 'l' && t[2] == 'm') return -10;
//     return 0;
// }


// int cost(UnorderedTriple<char> t) {
//     // // pyramid example below (3.1 + 3.11 are not sufficient) (3.4 is sufficient)
//     if (t[0] == 'a' && t[1] == 'b' && t[2] == 'e') return -75;
//     // pyramid example below (3.1 + 3.11 + 3.4 are not sufficient) (3.6 is sufficient for 10) (3.5 is sufficient for 100)
//     if (t[0] == 'b' && t[1] == 'c' && t[2] == 'd') return 10; // 10
//     if (t[0] == 'a' && t[1] == 'b' && t[2] == 'c') return -50;
//     if (t[0] == 'a' && t[1] == 'b' && t[2] == 'd') return -50;
//     if (t[0] == 'a' && t[1] == 'c' && t[2] == 'd') return -50;
//     return 0;
// }

int main() {
    std::vector<char> samples = {'a', 'b', 'c', 'd', 'e'};
    CubicSetPartitionProblem<char> problem(samples, cost);
    problem.solve();
    problem.printResults();
    
    return 0;
}