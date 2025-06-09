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
        Space::Vector r1 = u.generateOrthogonalVector().getNormalizedVector();
        Space::Vector r2 = u.crossProduct(r1).getNormalizedVector();
        std::cout << u.isOrthogonal(r1) << u.isOrthogonal(r2);
        std::cout << r1.isOrthogonal(u) << r1.isOrthogonal(r2);
        std::cout << r2.isOrthogonal(u) << r2.isOrthogonal(r1);
        std::cout << std::endl;
    }
    return 0;
}