#include <iostream>

#include "clustering_problem.hpp"
#include "instances.hpp"


int main() {
    ClusteringInstance<Space::Point> spaceInstance = generateSpaceInstance(2, 10, 100, 1);
    ClusteringProblem<Space::Point> problem(
        spaceInstance.unlabeledSamples,
        spaceInstance.cost
    );
    problem.solve();
    spaceInstance.printLabelEvaluation(problem.getLabels());
    return 0;
}