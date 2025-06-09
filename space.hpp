#ifndef SPACE_HPP
#define SPACE_HPP

#include <iostream>
#include <random>
#include <cfloat>
#include <cmath>


namespace Space {

    struct Vector {
        double x, y, z;

        Vector(double x, double y, double z);

        Vector operator*(double k) const;

        double operator*(Vector &other) const;

        bool isOrthogonal(Vector &other) const;

        static Vector generateUnitVector();
    };

    
    
}



#endif 