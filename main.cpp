#include <iostream>

#include "clustering_problem.hpp"
#include "instances.hpp"


int main() {
    const ClusteringInstance<Space::Point> SPACE_INSTANCE = generateSpaceInstance(2, 20, 100, 1);
    ClusteringProblem<Space::Point> problem(
        SPACE_INSTANCE.unlabeledSamples,
        SPACE_INSTANCE.cost
    );
    problem.solve();
    return 0;
}