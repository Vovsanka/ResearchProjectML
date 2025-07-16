#include "instances.hpp"

// define the samples
const std::vector<std::pair<char,int64_t>> SIMPLE_SAMPLES = {{'a',0}, {'b',1}, {'c',1}, {'d',1}};
const std::vector<std::pair<char,int64_t>> MULTICLUSTER_SAMPLES = {
    {'a',0}, {'b',0}, {'c',0}, {'d',0},
    {'e',1}, {'f',1}, {'g',1}, {'h',1}, {'i',1}
};
const std::vector<std::pair<char,int64_t>> PYRAMID_SAMPLES1 = {{'a',0}, {'b',0}, {'c',0}, {'d',0}, {'e',0}};
const std::vector<std::pair<char,int64_t>> PYRAMID_SAMPLES2 = {{'a',0}, {'b',0}, {'c',0}, {'d',0}};
const std::vector<std::pair<char,int64_t>> PYRAMID_SAMPLES_UNSOLVABLE = {{'a',0}, {'b',0}, {'c',0}, {'d',0}};

// predefine the cost functions
int64_t simpleCost(Utuple<3,char> t);
int64_t multiclusterCost(Utuple<3,char> t);
int64_t pyramidCost1(Utuple<3,char> t);
int64_t pyramidCost2(Utuple<3,char> t); 
int64_t pyramidCostUnsolvable(Utuple<3,char> t);

// define the clustering instances
const ClusteringInstance<char> SIMPLE_INSTANCE(SIMPLE_SAMPLES, simpleCost);
const ClusteringInstance<char> MULTICLUSTER_INSTANCE(MULTICLUSTER_SAMPLES, multiclusterCost);
const ClusteringInstance<char> PYRAMID_INSTANCE1(PYRAMID_SAMPLES1, pyramidCost1);
const ClusteringInstance<char> PYRAMID_INSTANCE2(PYRAMID_SAMPLES2, pyramidCost2);
const ClusteringInstance<char> PYRAMID_INSTANCE_UNSOLVABLE(PYRAMID_SAMPLES_UNSOLVABLE, pyramidCostUnsolvable);


ClusteringInstance<Space::Point> generateSpaceInstance(
    int64_t planeCount,
    int64_t pointsPerPlane,
    double maxDistance,
    double noise,
    unsigned int seed
) {
    std::vector<std::pair<Space::Point,int64_t>> labeledSamples = Space::generateSamplePointsOnDistinctPlanes(
        planeCount,
        pointsPerPlane,
        maxDistance,
        noise,
        seed
    );
    std::vector<Space::Point> points;
    for (auto [point, label] : labeledSamples) {
        points.push_back(point);
    }
    return ClusteringInstance<Space::Point>(
        labeledSamples,
        createSpaceCostFunction(points, maxDistance, noise)
    );
}

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
    if (t[0] == 'd' && t[1] == 'f' && t[2] == 'g') return -2;
    if (t[0] == 'f' && t[1] == 'g' && t[2] == 'i') return -10;
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

bool triangleIsLineLike(double a, double b, double c) {
    // compute the triangle angles [0, PI]
    double alpha = std::acos((b*b + c*c - a*a)/(2*b*c));
    double beta = std::acos((a*a + c*c - b*b)/(2*a*c));
    double gamma = std::acos((a*a + b*b - c*c)/(2*a*b));
    // inspect the largest angle
    double largestAngle = std::max(alpha, std::max(beta, gamma));
    return largestAngle > (120/180.0)*M_PI;
}

Space::Vector computeBestFittingPlaneNormalVector(std::vector<Space::Vector> locationVectors) {
    // the fitting plane contains the origin (no need to compute centroid)
    const int64_t N = locationVectors.size();
    // compute the scatter matrix
    std::array<std::array<double,3>,3> cov;
    for (int64_t i = 0; i < 3; i++) {
        for (int64_t j = 0; j < 3; j++) {
            cov[i][j] = 0;
            for (int64_t k = 0; k < N; k++) {
                cov[i][j] += locationVectors[k][i] * locationVectors[k][j];
            }
        }
    }
    // compute the eigenvector for the smallest eigenvalue (normalized norm of the best fitting plane)
    Eigen::Matrix3d A;
    A << cov[0][0], cov[0][1], cov[0][2],
         cov[1][0], cov[1][1], cov[1][2],
         cov[2][0], cov[2][1], cov[2][2];
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(A);
    if (solver.info() == Eigen::Success) {
        Eigen::Vector3d eigenvalues = solver.eigenvalues();
        Eigen::Matrix3d eigenvectors = solver.eigenvectors();
        int64_t minIndex;
        eigenvalues.minCoeff(&minIndex);
        Eigen::Vector3d smallestEigenvector = eigenvectors.col(minIndex);
        return Space::Vector(
            smallestEigenvector[0],
            smallestEigenvector[1],
            smallestEigenvector[2]
        ).getNormalizedVector();
    } 
    throw std::runtime_error("Eigen decomposition failed! ");
}

std::function<int64_t(Utuple<3,Space::Point>)> createSpaceCostFunction(
    const std::vector<Space::Point> &points,
    double maxDistance,
    double noise
) {
    const double TOL = maxDistance/1e4;
    std::function<double(Utuple<3,Space::Point>)> doubleCost = [TOL, points, maxDistance, noise](Utuple<3,Space::Point> pointTriple) -> double {
        // 0: compute triangle vectors and sides
        const int64_t pointCount = points.size();
        // compute the location vectors
        Space::Vector oa(pointTriple[0]);
        Space::Vector ob(pointTriple[1]);
        Space::Vector oc(pointTriple[2]);
        // compute the vectors of the triangle sides
        Space::Vector ab = ob - oa;
        Space::Vector bc = oc - ob;
        Space::Vector ca = oa - oc;
        // sort the sides
        std::array<double,3> sides = {ab.getLength(), bc.getLength(), ca.getLength()};
        std::sort(std::begin(sides), std::end(sides));
        // 1: skip if 2 triangle points are too close to each other (because of noise sensitivity)
        if (sides[0] < 0.5*maxDistance) return 0;
        // 2: skip line like triangles
        if (triangleIsLineLike(sides[0], sides[1], sides[2])) return 0;
        // 3: assign a penalty if the triangle is likely not from the same plane
        Space::Vector nBest = computeBestFittingPlaneNormalVector({oa, ob, oc});
        double ha = std::fabs(oa*nBest);
        double hb = std::fabs(ob*nBest);
        double hc = std::fabs(oc*nBest);
        if (ha + hb + hc > 3*noise + TOL) {
            return ((ha + hb + hc) - (3*noise + TOL))/(3*maxDistance);
        }
        // 4: skip triangles with too much noise because they span an unclear plane
        Space::Vector nTriangle = ab.crossProduct(ca*(-1)).getNormalizedVector();
        double ho = std::fabs(oa*nTriangle);
        if (ho > 10.0/pointCount*noise + TOL) return 0;
        // 5: compute the points that are likely in the same plane as the triangle points
        std::vector<Space::Vector> samePlaneVectors;
        for (auto &p : points) {
            double hp = std::fabs(Space::Vector(p)*nBest);
            if (hp < noise + TOL && Space::Vector(p).getLength() > 0.3*maxDistance) {
                samePlaneVectors.push_back(Space::Vector(p));
            }
        }
        Space::Vector nPlane = computeBestFittingPlaneNormalVector(samePlaneVectors);
        int64_t sameCount = 0;
        double reward = 0;
        for (auto &ov : samePlaneVectors) {
            double hv = std::fabs(ov*nPlane);
            double delta = (hv - (noise + TOL))/maxDistance;
            if (delta < 0) {
                reward += delta;
                sameCount++;
            }
        }
        if (sameCount <= 3) return 0;
        reward *= std::pow(2.0, sameCount - 4.0);
        return reward;
    };
    // double to int adapter
    return [TOL, doubleCost](Utuple<3,Space::Point> pointTriple) -> int64_t {
        double unroundedC = doubleCost(pointTriple);
        int64_t c = std::round(unroundedC/TOL);
        if (!c) return 0;
        return c;
    };
}
