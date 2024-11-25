#include "ray.h"

Ray::Ray() {}

Ray::Ray(const vec3& origin_, const vec3& direction_)
    : origin(origin_), direction(direction_) {}

vec3 Ray::at(double t) const {
    return origin + direction * t;
}