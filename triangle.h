#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "vec3.h"
#include "ray.h"
#include "material.h"

class Triangle {
public:
    vec3 v0, v1, v2;
    vec3 normal;
    Material* material_ptr;

    Triangle();
    Triangle(const vec3& v0_, const vec3& v1_, const vec3& v2_, Material* material_);

    bool rayIntersect(const Ray& ray, double t_min, double t_max, HitRecord& rec) const;
};

#endif

