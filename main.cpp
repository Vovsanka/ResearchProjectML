#include <iostream>
#include <set>
#include <map>
#include <vector>
#include <array>
#include <queue>
#include <functional>
#include <stdexcept>
#include <algorithm>
#include <boost/graph/push_relabel_max_flow.hpp>
#include <boost/graph/adjacency_list.hpp>


template <typename S = int> // domain S of samples
class UnorderedTriple { // unordered triple (sorted in the ascending order)
    std::array<S, 3> s;

public:

    UnorderedTriple(S s1, S s2, S s3) {
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
};


template <typename S = int> // domain S of samples
class CubicSetPartitionProblem {

    std::vector<S> samples; 
    std::function<int(UnorderedTriple<S>)> cost;  // cost function for triples SxSxS
    
    std::map<UnorderedTriple<>, int> tripleCosts;
    std::vector<std::vector<std::pair<int, int>>> relevantTriples;
    std::vector<std::vector<std::pair<int, int>>> negativeTriples;

    std::vector<int> indexClusterMapping;
 
    void initTripleCosts() {
        for (int i = 0; i < samples.size(); i++) {
            for (int j = i + 1; j < samples.size(); j++) {
                for (int k = j + 1; k < samples.size(); k++) {
                    UnorderedTriple<> indexTriple(i, j, k);
                    int c = getCost(getSampleTriple(indexTriple));
                    if (c) {
                        tripleCosts[indexTriple] = c;
                        relevantTriples[i].push_back(std::make_pair(j, k));
                        relevantTriples[j].push_back(std::make_pair(i, k));
                        relevantTriples[k].push_back(std::make_pair(i, j));
                    }
                    if (c < 0) {
                        negativeTriples[i].push_back(std::make_pair(j, k));
                        negativeTriples[j].push_back(std::make_pair(i, k));
                        negativeTriples[k].push_back(std::make_pair(i, j));
                    }
                }
            }
        }
    }

    UnorderedTriple<S> getSampleTriple(UnorderedTriple<> t) {
        return UnorderedTriple<S>(samples[t[0]], samples[t[1]], samples[t[2]]);
    }
 
    std::vector<CubicSetPartitionProblem<S>> getSubproblemsWithRegionGrowing() {
        std::vector<CubicSetPartitionProblem<S>> result;
        std::vector<bool> processed(samples.size(), false);
        for (int i = 0; i < samples.size(); i++) {
            if (processed[i]) continue;
            // start the BFS from i
            std::vector<int> chosen = {i};
            std::queue <int> q;
            processed[i] = true;
            q.push(i);
            while (!q.empty()) {
                int current = q.front();
                q.pop();
                for (auto [j, k] : negativeTriples[current]) {
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
            std::vector<S> subproblemSamples;
            for (auto idx : chosen) {
                subproblemSamples.push_back(samples[idx]);
            }
            result.push_back(CubicSetPartitionProblem<S>(subproblemSamples, cost));
        }
        return result;
    }
    
    
    
public: 
    explicit CubicSetPartitionProblem(const std::vector<S>& givenSamples, const std::function<int(UnorderedTriple<S>)> costCB)
    : samples(givenSamples), cost(costCB) {
        std::sort(std::begin(samples), std::end(samples)); // sort the samples in the ascending order
        relevantTriples.resize(samples.size(), {});
        negativeTriples.resize(samples.size(), {});
        indexClusterMapping.resize(samples.size(), 0);
        initTripleCosts();
    }

    int getCost(UnorderedTriple<S> t) {
        return cost(t);
    }

    int getCost(S s1, S s2, S s3) {
        return cost(UnorderedTriple<S>(s1, s2, s3));
    }

    std::map<S, int> getClusterMapping() {
        if (indexClusterMapping.empty()) throw std::runtime_error("The problem has not been solved yet!");
        std::map<S, int> result;
        for (int i = 0; i < samples.size(); i++) {
            result[samples[i]] = indexClusterMapping[i];
        }
        return result;
    } 

    void solve() {
        // init index of 
        std::map<S, int> indexOf;
        for (int i = 0; i < samples.size(); i++) {
            indexOf[samples[i]] = i;
        }
        std::vector<CubicSetPartitionProblem<S>> subproblems = getSubproblemsWithRegionGrowing();
        int clusterOffset = 0;
        for (auto& subproblem : subproblems) { // unite the solutions of the independent subproblems
            subproblem.solveAsSubproblem();
            int subclusterCount = 1;
            for (auto [sample, subcluster] : subproblem.getClusterMapping()) {
                indexClusterMapping[indexOf[sample]] = subcluster + clusterOffset;
                subclusterCount = std::max(subclusterCount, subcluster + 1);
            }
            clusterOffset += subclusterCount;
        }
    }

    void solveAsSubproblem() {
        // TODO !!!
        // just for now: all samples beint to the same cluster
    }
};


std::pair<int, std::vector<bool>> solveMinCut(int vertices, std::vector<std::tuple<int,int,int>> edges) {
    using namespace boost;

    // Define the graph type
    typedef adjacency_list<vecS, vecS, directedS,
        property<vertex_name_t, std::string>,
        property<edge_capacity_t, int,
            property<edge_residual_capacity_t, int,
                property<edge_reverse_t, adjacency_list<>::edge_descriptor>>>> Graph;

    int s = 0, t = vertices - 1;
            
    Graph g(vertices); 
    auto capacity = get(edge_capacity, g);
    auto rev = get(edge_reverse, g);
    auto residual_capacity = get(edge_residual_capacity, g);

    // Add edges with capacities
    for (auto& e : edges) {
        auto [u, v, c] = e;
        if (c < 0) throw std::runtime_error("MinCutProblem: negative edges are not allowed!");
        auto e1 = add_edge(u, v, g).first;
        auto e2 = add_edge(v, u, g).first; // Reverse edge
        capacity[e1] = capacity[e2] = c; // Reverse edge has 0 capacity
        rev[e1] = e2;
        rev[e2] = e1;
    };

    // Compute max flow
    int maxFlow = push_relabel_max_flow(g, s, t);

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

    return std::make_pair(maxFlow, visited);
}


int cost(UnorderedTriple<char> t) {
    if (t[0] == 'a' && t[1] == 'b' && t[2] == 'c')
        return -5;
    // if (t[0] == 'b' && t[1] == 'c' && t[2] == 'd')
    //     return -7;
    if (t[0] == 'a' && t[1] == 'b' && t[2] == 'd')
        return 10;
    return 0;
}


int main() {
    // std::vector<char> samples = {'a', 'b', 'c', 'd'};
    // CubicSetPartitionProblem<char> problem(samples, cost);
    // problem.solve();
    // for (auto [sample, cluster] : problem.getClusterMapping()) {
    //     std::cout << sample << " -> " << cluster << std::endl;
    // }
    //
    std::vector<std::tuple<int,int,int>> edges = {
        std::make_tuple(0, 1, 4),
        std::make_tuple(0, 2, 5),
        std::make_tuple(1, 3, 4),
        std::make_tuple(2, 3, 5),
        std::make_tuple(0, 3, 9) 
    };
    auto [maxFlow, partition] = solveMinCut(4, edges);
    std::cout << maxFlow << std::endl;
    for (int i = 0; i < 4; i++) {
        if (partition[i]) {
            std::cout << i << ' ';
        }
    }
    std::cout << std::endl;
    return 0;
}