#include <iostream>
#include <set>
#include <map>
#include <vector>
#include <optional>
#include <stdexcept>
#include <algorithm>


template <typename S = int> // domain S of samples
class CubicSetPartitionProblem {
    double (*cost)(S, S, S);  // cost function for triples SxSxS
    std::set<S> samples; 
    std::map<S, int> clusterMapping;

    std::optional<std::tuple<S,S,S>> getNegativeCostTuple(std::set<S> &subset) {
        std::set<S> rest;
        std::set_difference(
            samples.begin(), samples.end(),
            subset.begin(), subset.end(),
            std::inserter(rest, rest.begin())
        );
        for (S p : subset)
            for (S q : rest)
                for (S r : samples)
                    if (getCost(p, q, r) < 0)
                        return std::make_tuple(p, q, r);
        return std::nullopt;
    }

public: 
    explicit CubicSetPartitionProblem(std::set<S> samples, double (*costCB)(S, S, S)) {
        this->samples = samples;
        this->cost = costCB;
    }

    double getCost(S s1, S s2, S s3) {
        return cost(s1, s2, s3);
    }

    std::map<S, int> getClusterMapping() {
        if (clusterMapping.empty())
            throw std::runtime_error("The problem has not been solved!");
        return clusterMapping;
    } 

    int getCluster(S sample) {
        if (clusterMapping.empty())
            throw std::runtime_error("The problem has not been solved!");
        return clusterMapping[sample];
    }

    void solveWithRegionGrowing() {
        std::set<std::set<S>> R;
        std::set<S> Q;
        for (S sample : samples) {
            Q.insert(sample);
        }
        while (!Q.empty()) {
            auto popIterator = Q.begin();
            S p = *popIterator;
            Q.erase(popIterator);
            std::set<S> R1;
            R1.insert(p);
            do {   
                // TODO: optional function, pass R by reference!
                auto t = getNegativeCostTuple(R1); 
                if (t.has_value()) {
                    auto [p, q, r] = t.value();
                    R1.insert({p, q, r});
                    Q.erase(p);
                    Q.erase(q);
                    Q.erase(r);
                } else {
                    break;
                }
            } while (true);
            R.insert(R1);
        }
        // overwrite clusterMapping
        int clusterNumber = 0;
        for (std::set<S> cluster : R) {
            for (S sample : cluster) {
               clusterMapping[sample] = clusterNumber; 
            }
            clusterNumber++;
        }
    }
};


double cost(char s1, char s2, char s3) {
    std::vector<char> s = {s1, s2, s3};
    std::sort(s.begin(), s.end());
    if (s[0] == 'a' && s[1] == 'b' && s[2] == 'c')
        return -5;
    if (s[0] == 'b' && s[1] == 'c' && s[2] == 'd')
        return -7;
    if (s[0] == 'a' && s[1] == 'b' && s[2] == 'd')
        return 10;
    return 0;
}

int main() {
    std::set<char> samples = {'a', 'b', 'c', 'd'};
    CubicSetPartitionProblem<char> problem(samples, cost);
    problem.solveWithRegionGrowing();
    for (auto sample : samples) {
        std::cout << sample << " -> " << problem.getCluster(sample) << std::endl;
    }
    return 0;
}