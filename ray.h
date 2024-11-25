#ifndef RAY_H
#define RAY_H

#include "vec3.h"

class Ray {
public:
    vec3 origin;
    vec3 direction;

    Ray();
    Ray(const vec3& origin_, const vec3& direction_);

    vec3 at(double t) const;
};

#endif

