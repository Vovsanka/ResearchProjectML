#include <iostream>

#include "clustering_problem.hpp"
#include "instances.hpp"


int main() {
    ClusteringProblem<char> problem(
        MULTICLUSTER_INSTANCE.samples,
        MULTICLUSTER_INSTANCE.cost
    );
    problem.solve();
    return 0;
}