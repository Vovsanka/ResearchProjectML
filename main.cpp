#include <iostream>

#include "clustering_problem.hpp"
#include "instances.hpp"


int main() {
    // ClusteringProblem<char> problem(
    //     MULTICLUSTER_INSTANCE.samples,
    //     MULTICLUSTER_INSTANCE.cost
    // );
    // problem.solve();

    std::vector<Space::Plane> planes = Space::generateDistinctPlanes(3);
    for (auto &p : planes) {
        std::cout << p << std::endl;
    }

    return 0;
}