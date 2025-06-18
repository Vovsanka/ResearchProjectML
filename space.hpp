#ifndef SPACE_HPP
#define SPACE_HPP

#include <iostream>
#include <random>
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <string>
#include <sstream>
#include <iomanip>


namespace Space {

    const double TOL = 1e-6;

    struct Point {
        double x, y, z;
        std::string name;

        Point() = default;

        Point(double x, double y, double z, int64_t num);

        bool operator==(const Point &other) const;

        bool operator<(const Point &other) const;
    };

    struct Vector {
        double x, y, z;

        Vector() = default;

        Vector(double x, double y, double z);

        Vector(const Point &p);

        double operator[](int64_t i) const;

        bool operator==(const Vector &other) const;

        Vector operator+(const Vector &other) const;

        Vector operator-(const Vector &other) const;

        Vector operator*(double k) const;

        double operator*(const Vector &other) const;

        Vector crossProduct(const Vector &other) const;

        bool isOrthogonal(const Vector &other) const;

        bool isParallel(const Vector &other) const;

        double getLength() const;

        double getAngle(const Vector &other) const;

        Vector getNormalizedVector() const;

        static Vector generateUnitVector();

        Vector generateOrthogonalVector() const;
    };

    struct Plane {
        Vector n, r1, r2;

        Plane() = default;
        
        Plane(Vector norm);
        
        std::vector<Point> generatePoints(int64_t pointCount, int64_t startNum, double maxDistance, double maxNoise);
    };

    std::vector<Plane> generateDistinctPlanes(int64_t planeCount);

    std::vector<Point> generateSamplePointsOnDistinctPlanes(
        int64_t planeCount,
        int64_t pointsPerPlane,
        double maxDistance,
        double maxNoise
    );

    std::ostream& operator<<(std::ostream& os, const Vector &v);

    std::ostream& operator<<(std::ostream& os, const Point &p);

    std::ostream& operator<<(std::ostream& os, const Space::Plane &p);

}

#endif 