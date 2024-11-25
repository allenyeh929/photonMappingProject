#ifndef HIT_RECORD_H
#define HIT_RECORD_H

#include "vec3.h"
#include "ray.h"

class Material; // «e¦V«Å§i

struct HitRecord {
    vec3 point;
    vec3 normal;
    double t;
    bool front_face;
    Material* material_ptr;

    inline void set_face_normal(const Ray& r, const vec3& outward_normal) {
        front_face = vec3::dot(r.direction, outward_normal) < 0;
        normal = front_face ? outward_normal : -outward_normal;
    }
};

#endif

