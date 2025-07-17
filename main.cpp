#include <iostream>

#include "clustering_problem.hpp"
#include "instances.hpp"


template<typename S>
void solveEvaluateInstance(const ClusteringInstance<S> &instance, std::ofstream &csvData, bool costEvaluation = false) {
    if (costEvaluation) instance.evaluateCosts();
    ClusteringProblem<S> problem(
        instance.unlabeledSamples,
        instance.cost
    );
    auto startTime = std::chrono::high_resolution_clock::now();
    problem.solve();
    auto endTime = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> solvingDuration = endTime - startTime;
    std::cout << "Solving duration: " << std::fixed << std::setprecision(3) << solvingDuration.count() / 1e3 << " s" << std::endl;
    instance.printLabelEvaluation(problem.getLabels());
    auto [partialOptimality, accuracy] = instance.evaluateLabels(problem.getLabels());
    csvData << partialOptimality << "," << accuracy << ",";
    csvData << std::fixed << std::setprecision(3) << solvingDuration.count() / 1e3 << std::endl;
    std::cout << "----------------------------------------------" << std::endl;
}


int main() {
    // Pyramid instance (char) example
    // solveEvaluateInstance(PYRAMID_INSTANCE2, true);
    // Generate and solve subspace instances
    std::ifstream seedFile("../seeds.txt");
    if (!seedFile) {
        std::cerr << "Error opening the seed file!" << std::endl;
        return 1;
    }
    std::vector<unsigned int> seeds; 
    unsigned int value;
    while (seedFile >> value) {
        seeds.push_back(value);
    }
    seedFile.close();
    //
    std::ofstream csvData("data.csv");
    csvData << "opt,acc,time" << std::endl; 
    for (int64_t i = 0; i < 7; i++) {
        ClusteringInstance<Space::Point> spaceInstance = generateSpaceInstance(3, 20, 100, 3, seeds[i]);
        solveEvaluateInstance(spaceInstance, csvData);
    }
    csvData.close();
    //
    return 0;
}