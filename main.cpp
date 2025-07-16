#include <iostream>

#include "clustering_problem.hpp"
#include "instances.hpp"


int main() {
    // auto instance = PYRAMID_INSTANCE2;
    // ClusteringProblem<char> problem(
    //     instance.unlabeledSamples,
    //     instance.cost
    // );
    // auto startTime = std::chrono::high_resolution_clock::now();
    // problem.solve();
    // auto endTime = std::chrono::high_resolution_clock::now();
    // std::chrono::duration<double, std::milli> solvingDuration = endTime - startTime;
    // std::cout << "Solving duration: " << std::fixed << std::setprecision(3) << solvingDuration.count() / 1e3 << " s" << std::endl;
    // instance.printLabelEvaluation(problem.getLabels());
    // std::cout << "----------------------------------------------" << std::endl;
    // 
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
    for (int64_t i = 0; i < 15; i++) {
        ClusteringInstance<Space::Point> spaceInstance = generateSpaceInstance(3, 7, 100, 1, seeds[i]);
        // spaceInstance.evaluateCosts();
        ClusteringProblem<Space::Point> problem(
            spaceInstance.unlabeledSamples,
            spaceInstance.cost
        );
        auto startTime = std::chrono::high_resolution_clock::now();
        problem.solve();
        auto endTime = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> solvingDuration = endTime - startTime;
        std::cout << "Solving duration: " << std::fixed << std::setprecision(3) << solvingDuration.count() / 1e3 << " s" << std::endl;
        spaceInstance.printLabelEvaluation(problem.getLabels());
        std::cout << "----------------------------------------------" << std::endl;
    }
    // 
    return 0;
}