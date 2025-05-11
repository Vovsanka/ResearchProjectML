#include <iostream>
#include <set>
#include <map>
#include <vector>
#include <array>
#include <optional>
#include <stdexcept>
#include <algorithm>


template <typename S = int> // domain S of samples
class SampleTriple { // unordered triple (sorted in the ascending order)
    std::array<char, 3> s;

public:

    SampleTriple(S s1, S s2, S s3) {
        s = {s1, s2, s3};
        std::sort(s.begin(), s.end());
    }

    S operator[](int index) { // only getter, not setter (S&)
        return s[index];
    }
};


template <typename S = int> // domain S of samples
class CubicSetPartitionProblem {
    double (*cost)(SampleTriple<S>);  // cost function for triples SxSxS
    std::set<S> samples; 
    std::map<S, int> clusterMapping;

    std::optional<SampleTriple<S>> getNegativeCostTuple(std::set<S> &subset) {
        std::set<S> rest;
        std::set_difference(
            samples.begin(), samples.end(),
            subset.begin(), subset.end(),
            std::inserter(rest, rest.begin())
        );
        for (S p : subset)
            for (S q : rest)
                for (S r : samples) {
                    SampleTriple<S> t(p, q, r);
                    if (getCost(t) < 0)
                        return t;
                }
        return std::nullopt;
    }

public: 
    explicit CubicSetPartitionProblem(std::set<S> samples, double (*costCB)(SampleTriple<S>)) {
        this->samples = samples;
        this->cost = costCB;
    }

    double getCost(SampleTriple<S> t) {
        return cost(t);
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

    void applyRegionGrowing() {
        // TODO: region growing is only the first step
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
                auto opt_t = getNegativeCostTuple(R1); 
                if (opt_t.has_value()) {
                    SampleTriple<S> t = opt_t.value();
                    R1.insert({t[0], t[1], t[2]});
                    Q.erase(t[0]);
                    Q.erase(t[1]);
                    Q.erase(t[2]);
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

    void solve() {
        applyRegionGrowing();
    }
};


double cost(SampleTriple<char> t) {
    if (t[0] == 'a' && t[1] == 'b' && t[2] == 'c')
        return -5;
    if (t[0] == 'b' && t[1] == 'c' && t[2] == 'd')
        return -7;
    if (t[0] == 'a' && t[1] == 'b' && t[2] == 'd')
        return 10;
    return 0;
}

int main() {
    std::set<char> samples = {'a', 'b', 'c', 'd'};
    CubicSetPartitionProblem<char> problem(samples, cost);
    problem.solve();
    for (auto sample : samples) {
        std::cout << sample << " -> " << problem.getCluster(sample) << std::endl;
    }
    return 0;
}