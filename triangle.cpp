#include "triangle.h"
#include <cmath>

Triangle::Triangle() : material_ptr(nullptr) {}

Triangle::Triangle(const vec3& v0_, const vec3& v1_, const vec3& v2_, Material* material_)
    : v0(v0_), v1(v1_), v2(v2_), material_ptr(material_) {
    normal = (v1 - v0).cross(v2 - v0).normalize();
}

bool Triangle::rayIntersect(const Ray& ray, double t_min, double t_max, HitRecord& rec) const {
    const double EPSILON = 1e-8;
    vec3 edge1 = v1 - v0;
    vec3 edge2 = v2 - v0;
    vec3 h = ray.direction.cross(edge2);
    double a = vec3::dot(edge1, h);

    // 判斷光線是否平行於三角形
    if (std::fabs(a) < EPSILON)
        return false;

    double f = 1.0 / a;
    vec3 s = ray.origin - v0;
    double u = f * vec3::dot(s, h);
    if (u < 0.0 || u > 1.0)
        return false;

    vec3 q = s.cross(edge1);
    double v = f * vec3::dot(ray.direction, q);
    if (v < 0.0 || u + v > 1.0)
        return false;

    double t = f * vec3::dot(edge2, q);

    // 確保 t 在有效範圍內
    if (t < t_min || t > t_max)
        return false;

    // 填充 HitRecord
    rec.t = t;
    rec.point = ray.origin + ray.direction * t;
    rec.set_face_normal(ray, normal);
    rec.material_ptr = material_ptr; // 設置材質指標
    return true;
}

