#include "space.hpp"

using namespace Space;

std::mt19937 gen(std::random_device{}());

Vector::Vector(double x, double y, double z): x(x), y(y), z(z) {}

Vector Vector::operator*(double k) const {
    return Vector(k*x, k*y, k*z);
}

double Vector::operator*(const Vector &other) const {
    return x*other.x + y*other.y + z*other.z;
}

Vector Vector::crossProduct(const Vector &other) const {
    return Vector(
        y*other.z - z*other.y,
        z*other.x - x*other.z,
        x*other.y - y*other.x
    );
}

bool Vector::isOrthogonal(const Vector &other) const {
    return std::fabs((*this) * other) <= TOL;
}

double Vector::getLength() const {
    return sqrt(x*x + y*y + z*z);
}

Vector Vector::getNormalizedVector() const {
    double len = getLength();
    return Vector(
        x/len,
        y/len,
        z/len
    );
}

Vector Vector::generateUnitVector() {
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

Vector Vector::generateOrthogonalVector() const {
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

