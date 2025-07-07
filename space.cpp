#include "space.hpp"


// std::mt19937 gen(std::random_device{}());
std::mt19937 gen(42);

Space::Vector::Vector(double x, double y, double z): x(x), y(y), z(z) {}

Space::Vector::Vector(const Point &p) {
    x = p.x;
    y = p.y;
    z = p.z;
}

double Space::Vector::operator[](int64_t i) const {
    if (i == 0) return x;
    if (i == 1) return y;
    if (i == 2) return z;
    throw std::runtime_error("Vector component index out of range! ");
}

bool Space::Vector::operator==(const Vector &other) const {
    return std::fabs(x - other.x) <= TOL && std::fabs(y - other.y) <= TOL && std::fabs(z - other.z) <= TOL;
}

Space::Vector Space::Vector::operator+(const Vector &other) const {
    return Vector(
        x + other.x,
        y + other.y,
        z + other.z
    );
}

Space::Vector Space::Vector::operator-(const Vector &other) const {
    return Vector(
        x - other.x,
        y - other.y,
        z - other.z
    );
}

Space::Vector Space::Vector::operator*(double k) const {
    return Vector(k*x, k*y, k*z);
}

double Space::Vector::operator*(const Vector &other) const {
    return x*other.x + y*other.y + z*other.z;
}

Space::Vector Space::Vector::crossProduct(const Vector &other) const {
    return Vector(
        y*other.z - z*other.y,
        z*other.x - x*other.z,
        x*other.y - y*other.x
    );
}

bool Space::Vector::isOrthogonal(const Vector &other) const {
    return std::fabs((*this) * other) <= TOL; // scalar product is 0
}

bool Space::Vector::isParallel(const Vector &other) const {
    Vector cross = (*this).crossProduct(other); // cross product is 0 
    return std::fabs(cross.x) <= TOL && std::fabs(cross.y) <= TOL && std::fabs(cross.z) <= TOL;
}

double Space::Vector::getLength() const {
    return sqrt(x*x + y*y + z*z);
}

double Space::Vector::getAngle(const Vector &other) const {
    return std::acos(((*this)*other)/(getLength()*other.getLength()));
}

Space::Vector Space::Vector::getNormalizedVector() const {
    double len = getLength();
    return Vector(
        x/len,
        y/len,
        z/len
    );
}

Space::Vector Space::Vector::generateUnitVector() {
    std::uniform_real_distribution<double> thetaDist(0.0, 2.0 * M_PI); // [0, 2*PI)
    std::uniform_real_distribution<double> phiDist(-M_PI_2, std::nextafter(1.0, DBL_MAX)); // [-PI/2, +PI/2]
    double theta = thetaDist(gen);
    double phi = phiDist(gen);
    return Vector(
        std::cos(theta) * std::cos(phi),
        std::sin(theta) * std::cos(phi),
        std::sin(phi)
    );
}

Space::Vector Space::Vector::generateOrthogonalVector() const {
    double xo, yo, zo; // x*xo + y*yo + z*zo = 0
    // set the coordinates of all no[i] to 1 if n[i] = 0;
    std::vector<std::pair<double, double*>> allCoordinates = {
        std::make_pair(x, &xo),
        std::make_pair(y, &yo),
        std::make_pair(z, &zo),
    }, unsetCoordinates; 
    for (auto [c, coPtr] : allCoordinates) {
        if (std::fabs(c) <= TOL) *coPtr = 1;
        else unsetCoordinates.push_back(std::make_pair(c, coPtr));
    }
    if (unsetCoordinates.empty()) std::runtime_error("Error generating an orthogonal vector. ");
    // set all unset coordinates except of the first one to 1, recompute the rhs
    double rhs = 0; 
    for (int64_t i = 1; i < unsetCoordinates.size(); i++) {
        auto [c, coPtr] = unsetCoordinates[i];
        *coPtr = 1;
        rhs -= c;
    }
    // compute the first unset coordinate sucht that the resulting vector is orthogonal
    auto [c, coPtr] = unsetCoordinates[0];
    *coPtr = rhs/c;
    return Vector(xo, yo, zo);
}

Space::Point::Point(double x, double y, double z, int64_t num) : x(x), y(y), z(z) {
    std::ostringstream os;
    int64_t letter = num%52; // 52 = 26*2
    if (letter < 26) {
        os << char('a' + letter);
    } else {
        os << char('A' + (letter - 26));
    }
    if (num >= 52) {
        os << num/52;
    }
    name = os.str();
}


bool Space::Point::operator==(const Point &other) const {
    return (name == other.name);
}

bool Space::Point::operator<(const Point &other) const {
    return (name < other.name);
}

Space::Plane::Plane(Vector norm) {
    n = norm.getNormalizedVector();
    r1 = n.generateOrthogonalVector().getNormalizedVector();
    r2 = n.crossProduct(r1).getNormalizedVector();
}

std::vector<Space::Point> Space::Plane::generatePoints(
    int64_t pointCount,
    int64_t startNum,
    double maxDistance,
    double noise
) {
    std::uniform_real_distribution<double> planeDist(-maxDistance, std::nextafter(maxDistance, DBL_MAX)); // [-maxDistance, +maxDistance]
    std::normal_distribution<double> noiseDist(0.0, noise);
    std::vector<Point> points;
    for (int64_t i = 0; i < pointCount; i++) {
        double k1 = planeDist(gen);
        double k2 = planeDist(gen);
        double kn = noiseDist(gen);
        Vector positionVector = r1*k1 + r2*k2 + n*kn;
        points.push_back(Point(
            positionVector.x,
            positionVector.y,
            positionVector.z,
            startNum + i
        ));
    }
    return points;
}

std::vector<Space::Plane> Space::generateDistinctPlanes(int64_t planeCount) {
    const double LOWER_ANGLE_ALMOST_PARALLEL = 45;
    const double UPPER_ANGLE_ALMOST_PARALLEL = 180 - LOWER_ANGLE_ALMOST_PARALLEL;
    std::vector<Vector> norms;
    for (int64_t i = 0; i < planeCount; i++) {
        bool denied;
        Vector candidate;
        do {
            denied = false;
            candidate = Vector::generateUnitVector();
            for (Vector &n : norms) {
                double angle = candidate.getAngle(n)/M_PI * 180.0;
                if (angle < LOWER_ANGLE_ALMOST_PARALLEL || UPPER_ANGLE_ALMOST_PARALLEL < angle) {
                    denied = true;
                    break;
                }
            }
        } while (denied);
        norms.push_back(candidate);
    }
    // build the planes
    std::vector<Plane> planes;
    for (auto &n : norms) {
        planes.push_back(Plane(n));
    }
    return planes;
}


std::vector<std::pair<Space::Point,int64_t>> Space::generateSamplePointsOnDistinctPlanes(
    int64_t planeCount,
    int64_t pointsPerPlane,
    double maxDistance,
    double noise
) {
    std::ofstream csvPoints("points.csv");
    std::ofstream csvPlanes("planes.csv");
    csvPoints << "p,x,y,z" << std::endl;
    csvPlanes << "p,x,y,z" << std::endl;
    std::vector<Space::Plane> planes = Space::generateDistinctPlanes(planeCount);
    std::cout << "Generating cubic space clustering instance: " << std::endl;
    std::vector<std::pair<Space::Point,int64_t>> samples;
    int64_t startNum = 0;
    for (int64_t i = 0; i < planeCount; i++) {
        std::cout << planes[i] << "\nPoints: ";
        std::vector<Space::Point> points = planes[i].generatePoints(pointsPerPlane, startNum, maxDistance, noise);
        csvPlanes << i << "," << planes[i].n.x << "," << planes[i].n.y << "," << planes[i].n.z << std::endl; 
        for (auto &p : points) {
            samples.push_back(std::make_pair(p, i));
            std::cout << p << " ";
            csvPoints << i << "," << p.x << "," << p.y << "," << p.z << std::endl;
        }
        std::cout << "\n" << std::endl;
        startNum += pointsPerPlane;
    }
    csvPoints.close();
    csvPlanes.close();
    std::shuffle(std::begin(samples), std::end(samples), gen); // random shuffle the samples
    return samples;
}

std::ostream& Space::operator<<(std::ostream& os, const Vector &v) {
    os << std::fixed << std::setprecision(2);
    os << "(" << v.x << ", " << v.y << ", " << v.z << ")";
    return os;
}

std::ostream& Space::operator<<(std::ostream& os, const Point &p) {
    os << p.name;
    return os;
}

std::ostream& Space::operator<<(std::ostream& os, const Space::Plane &p) {
    os << "Plane defined by:\n";
    os << "n=" << p.n << "\n";
    os << "r1=" << p.r1 << "\n";
    os << "r2=" << p.r2 << "\n";
    return os;
}