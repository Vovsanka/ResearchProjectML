#include <iostream>

#include "clustering_problem.hpp"
#include "instances.hpp"


int main() {
    const ClusteringInstance<Space::Point> SPACE_INSTANCE = generateSpaceInstance(3, 10, 100, 0);
    ClusteringProblem<Space::Point> problem(
        SPACE_INSTANCE.samples,
        SPACE_INSTANCE.cost
    );
    problem.solve();
    return 0;
}