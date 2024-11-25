#include "sphere.h"
#include <cmath>

Sphere::Sphere() {}

Sphere::Sphere(const vec3& center_, double radius_, Material* material_)
    : center(center_), radius(radius_), material_ptr(material_) {}

bool Sphere::rayIntersect(const Ray& ray, double t_min, double t_max, HitRecord& rec) const {
    vec3 oc = ray.origin - center;
    double a = ray.direction.length_square();
    double half_b = vec3::dot(oc, ray.direction);
    double c = oc.length_square() - radius * radius;
    double discriminant = half_b * half_b - a * c;

    if (discriminant > 0) {
        double sqrt_d = std::sqrt(discriminant);
        double root = (-half_b - sqrt_d) / a;
        if (root < t_max && root > t_min) {
            rec.t = root;
            rec.point = ray.at(rec.t);
            vec3 outward_normal = (rec.point - center) / radius;
            rec.set_face_normal(ray, outward_normal);
            rec.material_ptr = material_ptr;
            return true;
        }
        root = (-half_b + sqrt_d) / a;
        if (root < t_max && root > t_min) {
            rec.t = root;
            rec.point = ray.at(rec.t);
            vec3 outward_normal = (rec.point - center) / radius;
            rec.set_face_normal(ray, outward_normal);
            rec.material_ptr = material_ptr;
            return true;
        }
    }
    return false;
}
