#include <iostream>

#include "clustering_problem.hpp"
#include "instances.hpp"


int main() {
    std::vector<Space::Point> samples = {
            Space::Point(5, 5, 0, 0),
            Space::Point(1, 2, 0, 1),
            Space::Point(-5, -5, 0, 2),
            Space::Point(5, 4, 0, 3), 
            Space::Point(0, 0, 10, 4),
            Space::Point(0, 0, 5, 5),
            Space::Point(0, 0, -10, 6)
        };
    ClusteringProblem<Space::Point> problem(
        CUBIC_SPACE_INSTANCE.samples,
        // samples,
        CUBIC_SPACE_INSTANCE.cost
    );
    problem.solve();
    return 0;
}