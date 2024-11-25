#ifndef SPHERE_H
#define SPHERE_H

#include "vec3.h"
#include "ray.h"
#include "material.h"
#include "hit_record.h"

class Sphere {
public:
    vec3 center;
    double radius;
    Material* material_ptr;

    Sphere();
    Sphere(const vec3& center_, double radius_, Material* material_);

    bool rayIntersect(const Ray& ray, double t_min, double t_max, HitRecord& rec) const;
};

#endif

