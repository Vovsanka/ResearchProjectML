#ifndef INSTANCES_HPP
#define INSTANCES_HPP

#include <iostream>
#include <vector>
#include <array>
#include <map>
#include <functional>

#include <Eigen/Dense>

#include "clustering_problem.hpp"
#include "space.hpp"


template <typename S>
struct ClusteringInstance {
    std::vector<S> unlabeledSamples;
    std::vector<int64_t> actualClustering;
    std::function<int64_t(Utuple<3,S>)> cost;
    std::function<int64_t(Utuple<2,S>)> pairCost;

    ClusteringInstance(
        std::vector<std::pair<S,int64_t>> labeledSamples,
        std::function<int64_t(Utuple<3,S>)> tripleCostCB,
        std::function<int64_t(Utuple<2,S>)> pairCostCB= [](Utuple<2,S> p)->int64_t{return 0;}
    ) : cost(tripleCostCB), pairCost(pairCostCB) {
        this->actualClustering = std::vector<int64_t>(labeledSamples.size());
        for (int64_t i = 0; i < labeledSamples.size(); i++) {
            this->unlabeledSamples.push_back(labeledSamples[i].first);
            this->actualClustering[i] = labeledSamples[i].second;
        }
    }

    std::pair<double,double> evaluateLabels(const std::map<Upair,int64_t> &labels) const {
        // computes partial optimality percentage and accuracy
        int64_t sampleCount = unlabeledSamples.size();
        int64_t labelCount = sampleCount*(sampleCount - 1);
        int64_t fixedLabelCount = 0;
        int64_t tp = 0, tn = 0, fp = 0, fn = 0;
        for (int64_t i = 0; i < sampleCount; i++) {
            for (int64_t j = i + 1; j < sampleCount; j++) {
                labelCount++;
                int64_t l = labels.at(Upair({i, j}));
                if (l > 0) {
                    fixedLabelCount++;
                    if (l == 1) { 
                        // check cut correctness (negatives)
                        if (actualClustering[i] != actualClustering[j]) tn++;
                        else fn++;
                    }
                    if (l == 2) {
                        // check join correctness (positives)
                        if (actualClustering[i] == actualClustering[j]) tp++;
                        else fp++;
                    }
                }
            }
        }
        double partialOptimality = (1.0 * fixedLabelCount) / labelCount;
        double accuracy = (1.0 * (tp + tn)) / (tp + tn + fp + fn);
        return std::make_pair(partialOptimality, accuracy);
    }

    void printLabelEvaluation(const std::map<Upair,int64_t> &labels) const {
        auto [partialOptimality, accuracy] = evaluateLabels(labels);
        std::cout << "Partial optimality (% of the fixed labels): ";
        std::cout << std::fixed << std::setprecision(1) << partialOptimality * 100 << "%\n";
        std::cout << "Accuracy (for the fixed labels): ";
        std::cout << std::fixed << std::setprecision(1) << accuracy * 100 << "%";
        std::cout << std::endl;
    }

    void evaluateCosts() {
        std::vector<int64_t> same, diff;
        int64_t sampleCount = unlabeledSamples.size();
        for (int64_t i = 0; i < sampleCount; i++) {
            for (int64_t j = i + 1; j < sampleCount; j++) {
                int64_t c = pairCost(Utuple<2,S>({unlabeledSamples[i], unlabeledSamples[j]}));
                if (actualClustering[i] == actualClustering[j]) {
                    same.push_back(c);
                } else {
                    diff.push_back(c);
                }
            }
        }
        for (int64_t i = 0; i < sampleCount; i++) {
            for (int64_t j = i + 1; j < sampleCount; j++) {
                for (int64_t k = j + 1; k < sampleCount; k++) {
                    int64_t c = cost(Utuple<3,S>({unlabeledSamples[i], unlabeledSamples[j], unlabeledSamples[k]}));
                    if (!c) continue;
                    if (actualClustering[i] == actualClustering[j] && actualClustering[i] == actualClustering[k]) {
                        same.push_back(c);
                    } else {
                        diff.push_back(c);
                    }
                }
            }
        }
        std::sort(std::begin(same), std::end(same));
        std::sort(std::begin(diff), std::end(diff));
        std::ofstream sameFile("same.txt"), diffFile("diff.txt");
        for (auto c : same) {
            sameFile << c << std::endl;
        }
        for (auto c : diff) {
            diffFile << c << std::endl;
        }
        sameFile.close();
        diffFile.close();
    }

};


ClusteringInstance<Space::Point> generateSpaceInstance(
    int64_t planeCount,
    int64_t pointsPerPlane,
    double maxDistance,
    double maxNoise
);

std::function<int64_t(Utuple<3,Space::Point>)> createSpaceCostFunction(
    const std::vector<Space::Point> &points,
    double maxDistance,
    double maxNoise
);

extern const ClusteringInstance<char> SIMPLE_INSTANCE;
extern const ClusteringInstance<char> MULTICLUSTER_INSTANCE;
extern const ClusteringInstance<char> PYRAMID_INSTANCE1;
extern const ClusteringInstance<char> PYRAMID_INSTANCE2; 
extern const ClusteringInstance<char> PYRAMID_INSTANCE_UNSOLVABLE;

#endif