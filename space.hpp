#ifndef SPACE_HPP
#define SPACE_HPP

#include <iostream>
#include <random>
#include <cfloat>
#include <cmath>


namespace Space {

    const double TOL = 1e-6;

    struct Vector {
        double x, y, z;

        Vector() = default;

        Vector(double x, double y, double z);

        Vector operator*(double k) const;

        double operator*(const Vector &other) const;

        Vector crossProduct(const Vector &other) const;

        bool isOrthogonal(const Vector &other) const;

        double getLength() const;

        Vector getNormalizedVector() const;

        static Vector generateUnitVector();

        Vector generateOrthogonalVector() const;
    };

    struct Plane {
        Vector n, r1, r2;
        
        Plane(Vector norm);
        
        // Point generatePoints(int64_t pointCount, noise)
    };

    std::vector<Plane> generateDistinctPlanes(int64_t planeCount = 3);

    std::ostream& operator<<(std::ostream& os, const Vector &v);

    std::ostream& operator<<(std::ostream& os, const Space::Plane &p);

}

#endif 