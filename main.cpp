#include <iostream>
#include <set>
#include <map>
#include <vector>
#include <array>
#include <queue>
#include <functional>
#include <stdexcept>
#include <algorithm>


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
    std::function<double(UnorderedTriple<S>)> cost;  // cost function for triples SxSxS
    
    std::map<UnorderedTriple<>, double> tripleCosts;
    std::vector<std::vector<std::pair<int, int>>> relevantTriples;
    std::vector<std::vector<std::pair<int, int>>> negativeTriples;

    std::vector<int> indexClusterMapping;
 
    void initTripleCosts() {
        for (int i = 0; i < samples.size(); i++) {
            for (int j = i + 1; j < samples.size(); j++) {
                for (int k = j + 1; k < samples.size(); k++) {
                    UnorderedTriple<> indexTriple(i, j, k);
                    double c = getCost(getSampleTriple(indexTriple));
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
    explicit CubicSetPartitionProblem(const std::vector<S>& givenSamples, const std::function<double(UnorderedTriple<S>)> costCB)
    : samples(givenSamples), cost(costCB) {
        std::sort(std::begin(samples), std::end(samples)); // sort the samples in the ascending order
        relevantTriples.resize(samples.size(), {});
        negativeTriples.resize(samples.size(), {});
        indexClusterMapping.resize(samples.size(), 0);
        initTripleCosts();
    }

    double getCost(UnorderedTriple<S> t) {
        return cost(t);
    }

    double getCost(S s1, S s2, S s3) {
        return cost(UnorderedTriple<S>(s1, s2, s3));
    }

    std::map<S, int> getClusterMapping() {
        if (indexClusterMapping.empty())
            throw std::runtime_error("The problem has not been solved yet!");
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
        // just for now: all samples belong to the same cluster
    }
};


double cost(UnorderedTriple<char> t) {
    if (t[0] == 'a' && t[1] == 'b' && t[2] == 'c')
        return -5;
    // if (t[0] == 'b' && t[1] == 'c' && t[2] == 'd')
    //     return -7;
    if (t[0] == 'a' && t[1] == 'b' && t[2] == 'd')
        return 10;
    return 0;
}


int main() {
    std::vector<char> samples = {'a', 'b', 'c', 'd'};
    CubicSetPartitionProblem<char> problem(samples, cost);
    problem.solve();
    for (auto [sample, cluster] : problem.getClusterMapping()) {
        std::cout << sample << " -> " << cluster << std::endl;
    }
    return 0;
}