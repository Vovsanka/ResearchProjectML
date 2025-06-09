#include "instances.hpp"

// define the samples
const std::vector<char> SIMPLE_SAMPLES = {'a', 'b', 'c', 'd'};
const std::vector<char> MULTICLUSTER_SAMPLES = {'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm'};
const std::vector<char> PYRAMID_SAMPLES = {'a', 'b', 'c', 'd', 'e'};

// predefine the cost functions
int64_t simple_cost(Utuple<3,char> t);
int64_t multicluster_cost(Utuple<3,char> t);
int64_t pyramid_cost1(Utuple<3,char> t);
int64_t pyramid_cost2(Utuple<3,char> t); 
int64_t pyramid_cost_unsolvable(Utuple<3,char> t);

// define the clustering instances
const ClusteringInstance<char> SIMPLE_INSTANCE(SIMPLE_SAMPLES, simple_cost);
const ClusteringInstance<char> MULTICLUSTER_INSTANCE(MULTICLUSTER_SAMPLES, multicluster_cost);
const ClusteringInstance<char> PYRAMID_INSTANCE1(PYRAMID_SAMPLES, pyramid_cost1);
const ClusteringInstance<char> PYRAMID_INSTANCE2(PYRAMID_SAMPLES, pyramid_cost2);
const ClusteringInstance<char> PYRAMID_INSTANCE_UNSOLVABLE(PYRAMID_SAMPLES, pyramid_cost_unsolvable);

int64_t simple_cost(Utuple<3,char> t) {
    if (t[0] == 'a' && t[1] == 'b' && t[2] == 'c')
        return -5;
    if (t[0] == 'b' && t[1] == 'c' && t[2] == 'd')
        return -22;
    if (t[0] == 'a' && t[1] == 'b' && t[2] == 'd')
        return 10;
    return 0;
}

int64_t multicluster_cost(Utuple<3,char> t) {
    // example: 3.1 and 3.11 are sufficient
    if (t[0] == 'a' && t[1] == 'b' && t[2] == 'c') return -1;
    if (t[0] == 'a' && t[1] == 'c' && t[2] == 'd') return -15;
    if (t[0] == 'd' && t[1] == 'e' && t[2] == 'h') return 50;
    if (t[0] == 'e' && t[1] == 'f' && t[2] == 'h') return -50;
    if (t[0] == 'f' && t[1] == 'g' && t[2] == 'i') return -30;
    if (t[0] == 'd' && t[1] == 'h' && t[2] == 'k') return -2;
    if (t[0] == 'i' && t[1] == 'k' && t[2] == 'l') return -4;
    if (t[0] == 'j' && t[1] == 'l' && t[2] == 'm') return -10;
    return 0;
}

int64_t pyramid_cost1(Utuple<3,char> t) {
    // // pyramid example (3.1 + 3.11 are not sufficient) (3.4 is sufficient)
    if (t[0] == 'a' && t[1] == 'b' && t[2] == 'e') return -75;
    if (t[0] == 'b' && t[1] == 'c' && t[2] == 'd') return 10; 
    if (t[0] == 'a' && t[1] == 'b' && t[2] == 'c') return -50;
    if (t[0] == 'a' && t[1] == 'b' && t[2] == 'd') return -50;
    if (t[0] == 'a' && t[1] == 'c' && t[2] == 'd') return -50;
    return 0;
}

int64_t pyramid_cost2(Utuple<3,char> t) {
    // pyramid example below (3.1 + 3.11 + 3.4 are not sufficient) (3.6 is sufficient for 10)
    if (t[0] == 'b' && t[1] == 'c' && t[2] == 'd') return 10; 
    if (t[0] == 'a' && t[1] == 'b' && t[2] == 'c') return -50;
    if (t[0] == 'a' && t[1] == 'b' && t[2] == 'd') return -50;
    if (t[0] == 'a' && t[1] == 'c' && t[2] == 'd') return -50;
    return 0;
}

int64_t pyramid_cost_unsolvable(Utuple<3,char> t) {
    // // pyramid example (everything is insufficient, 3.3 separates the pyramid base for the cost 100+)
    if (t[0] == 'b' && t[1] == 'c' && t[2] == 'd') return 100;
    if (t[0] == 'a' && t[1] == 'b' && t[2] == 'c') return -50;
    if (t[0] == 'a' && t[1] == 'b' && t[2] == 'd') return -50;
    if (t[0] == 'a' && t[1] == 'c' && t[2] == 'd') return -50;
    return 0;
}


