#include "instances.hpp"

// define the samples
const std::vector<char> SIMPLE_SAMPLES = {'a', 'b', 'c', 'd'};
const std::vector<char> MULTICLUSTER_SAMPLES = {'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm'};
const std::vector<char> PYRAMID_SAMPLES = {'a', 'b', 'c', 'd', 'e'};

// predefine the cost functions
int64_t simpleCost(Utuple<3,char> t);
int64_t multiclusterCost(Utuple<3,char> t);
int64_t pyramidCost1(Utuple<3,char> t);
int64_t pyramidCost2(Utuple<3,char> t); 
int64_t pyramidCostUnsolvable(Utuple<3,char> t);
double cubicSpaceCost(Utuple<3,Space::Point> t);

// define the clustering instances
const ClusteringInstance<char> SIMPLE_INSTANCE(SIMPLE_SAMPLES, simpleCost);
const ClusteringInstance<char> MULTICLUSTER_INSTANCE(MULTICLUSTER_SAMPLES, multiclusterCost);
const ClusteringInstance<char> PYRAMID_INSTANCE1(PYRAMID_SAMPLES, pyramidCost1);
const ClusteringInstance<char> PYRAMID_INSTANCE2(PYRAMID_SAMPLES, pyramidCost2);
const ClusteringInstance<char> PYRAMID_INSTANCE_UNSOLVABLE(PYRAMID_SAMPLES, pyramidCostUnsolvable);
const ClusteringInstance<Space::Point> CUBIC_SPACE_INSTANCE(
    Space::generateSamplePointsOnDistinctPlanes(2, 30, 10, 0),
    doubleToIntCostWrapper<Utuple<3,Space::Point>>(cubicSpaceCost)
);

int64_t simpleCost(Utuple<3,char> t) {
    if (t[0] == 'a' && t[1] == 'b' && t[2] == 'c')
        return -5;
    if (t[0] == 'b' && t[1] == 'c' && t[2] == 'd')
        return -22;
    if (t[0] == 'a' && t[1] == 'b' && t[2] == 'd')
        return 10;
    return 0;
}

int64_t multiclusterCost(Utuple<3,char> t) {
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

int64_t pyramidCost1(Utuple<3,char> t) {
    // // pyramid example (3.1 + 3.11 are not sufficient) (3.4 is sufficient)
    if (t[0] == 'a' && t[1] == 'b' && t[2] == 'e') return -75;
    if (t[0] == 'b' && t[1] == 'c' && t[2] == 'd') return 10; 
    if (t[0] == 'a' && t[1] == 'b' && t[2] == 'c') return -50;
    if (t[0] == 'a' && t[1] == 'b' && t[2] == 'd') return -50;
    if (t[0] == 'a' && t[1] == 'c' && t[2] == 'd') return -50;
    return 0;
}

int64_t pyramidCost2(Utuple<3,char> t) {
    // pyramid example below (3.1 + 3.11 + 3.4 are not sufficient) (3.6 is sufficient for 10)
    if (t[0] == 'b' && t[1] == 'c' && t[2] == 'd') return 10; 
    if (t[0] == 'a' && t[1] == 'b' && t[2] == 'c') return -50;
    if (t[0] == 'a' && t[1] == 'b' && t[2] == 'd') return -50;
    if (t[0] == 'a' && t[1] == 'c' && t[2] == 'd') return -50;
    return 0;
}

int64_t pyramidCostUnsolvable(Utuple<3,char> t) {
    // // pyramid example (everything is insufficient, 3.3 separates the pyramid base for the cost 100+)
    if (t[0] == 'b' && t[1] == 'c' && t[2] == 'd') return 100;
    if (t[0] == 'a' && t[1] == 'b' && t[2] == 'c') return -50;
    if (t[0] == 'a' && t[1] == 'b' && t[2] == 'd') return -50;
    if (t[0] == 'a' && t[1] == 'c' && t[2] == 'd') return -50;
    return 0;
}

double cubicSpaceCost(Utuple<3,Space::Point> t) {
    // compute the triangle sides and perimeter
    double a = t[0].getDistance(t[1]);
    double b = t[1].getDistance(t[2]);
    double c = t[2].getDistance(t[0]);
    double p = a + b + c;
    // no cost if one 2 triangle points are too close to each other
    std::array<double, 3> sides = {a, b, c};
    std::sort(std::begin(sides), std::end(sides)); 
    if (sides[0] < (p/3)/2 ) return 0;
    // compute the triangle angles [0, PI]
    double alpha = std::acos((b*b + c*c - a*a)/(2*b*c));
    double beta = std::acos((a*a + c*c - b*b)/(2*a*c));
    double gamma = std::acos((a*a + b*b - c*c)/(2*a*b));
    // inspect the largest angle
    double largestAngle = std::max(alpha, std::max(beta, gamma));
    const double REWARD_BOUND = (175/180.0)*M_PI;
    const double PENALTY_BOUND = (65/180.0)*M_PI;
    if (largestAngle > REWARD_BOUND) return -(largestAngle - REWARD_BOUND)*(largestAngle - REWARD_BOUND) * p;
    if (largestAngle < PENALTY_BOUND) return sqrt((PENALTY_BOUND - largestAngle) * p);
    return 0;
}
