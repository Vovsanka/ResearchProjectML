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
    for (auto &plane : planes) {
        std::cout << plane;
        std::cout << plane.n.getLength() << plane.r1.getLength() << plane.r2.getLength();
        std::cout << plane.n.isOrthogonal(plane.r1);
        std::cout << plane.n.isOrthogonal(plane.r2);
        std::cout << plane.r1.isOrthogonal(plane.r2);
        std::cout << "\n" << std::endl;
    }

    return 0;
}