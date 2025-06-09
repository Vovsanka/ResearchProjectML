#include "space.hpp"

using namespace Space;

std::mt19937 gen(std::random_device{}());

Vector::Vector(double x, double y, double z): x(x), y(y), z(z) {}

Vector Vector::operator*(double k) const {
    return Vector(k*x, k*y, k*z);
}

double Vector::operator*(Vector &other) const {
    return x*other.x + y*other.y + z*other.z;
}

bool Vector::isOrthogonal(Vector &other) const {
    return abs((*this) * other) < 1e-6;
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




