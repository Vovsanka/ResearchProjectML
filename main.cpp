#include <iostream>

#include "clustering_problem.hpp"
#include "instances.hpp"


int main() {
    ClusteringProblem<Space::Point> problem(
        CUBIC_SPACE_INSTANCE.samples,
        CUBIC_SPACE_INSTANCE.cost
    );
    problem.solve();
    return 0;
}