#include <iostream>

#include "clustering_problem.hpp"
#include "instances.hpp"


int main() {
    // ClusteringProblem<char> problem(
    //     MULTICLUSTER_INSTANCE.samples,
    //     MULTICLUSTER_INSTANCE.cost
    // );
    // problem.solve();
    for (int i = 0; i < 10; i++) {
        Space::Vector u = Space::Vector::generateUnitVector();
        std::cout << sqrt(u.x*u.x + u.y*u.y + u.z*u.z) << ": ";
        std::cout << u.x << " " << u.y << " " << u.z << std::endl;
    }
    return 0;
}