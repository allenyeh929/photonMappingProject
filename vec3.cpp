#include "vec3.h"

#include <cmath>
#include <algorithm>

vec3::vec3() : x(0), y(0), z(0) {}

vec3::vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}

double& vec3::operator[](int index) {
    if (index == 0) return x;
    else if (index == 1) return y;
    else if (index == 2) return z;
}

const double& vec3::operator[](int index) const {
    if (index == 0) return x;
    else if (index == 1) return y;
    else if (index == 2) return z;
}

vec3 vec3::operator+(const vec3& v) const {
    return vec3(x + v.x, y + v.y, z + v.z);
}

vec3 vec3::operator-(const vec3& v) const {
    return vec3(x - v.x, y - v.y, z - v.z);
}

vec3 vec3::operator*(double s) const {
    return vec3(x * s, y * s, z * s);
}

vec3 vec3::operator/(double s) const {
    return vec3(x / s, y / s, z / s);
}

vec3 vec3::operator-() const {
    return vec3(-x, -y, -z);
}

// 向量與向量的逐元素相乘
vec3 vec3::operator*(const vec3& other) const {
    return vec3(x * other.x, y * other.y, z * other.z);
}

// 向量與純量的乘法賦值
vec3& vec3::operator*=(double scalar) {
    x *= scalar;
    y *= scalar;
    z *= scalar;
    return *this;
}

vec3 vec3::operator/(const vec3& other) const {
    return vec3(x / other.x, y / other.y, z / other.z);
}

double vec3::dot(const vec3 & u, const vec3 & v) {
    return u.x * v.x
        + u.y * v.y
        + u.z * v.z;
}

vec3 vec3::cross(const vec3& v) const {
    return vec3(
        y * v.z - z * v.y,
        z * v.x - x * v.z,
        x * v.y - y * v.x
    );
}

double vec3::length() const {
    return std::sqrt(length_square());
}

double vec3::length_square() const {
    return x * x + y * y + z * z;
}

bool vec3::near_zero() const {
    // Return true if the vector is close to zero in all dimensions.
    auto s = 1e-8;
    return (std::fabs(x) < s) && (std::fabs(y) < s) && (std::fabs(z) < s);
}

vec3 vec3::normalize() const {
    double len = length();
    if (len > 0)
        return (*this) / len;
    else
        return *this;
}

vec3 vec3::reflect(const vec3& I, const vec3& N) {
    return I - N * 2.0 * dot(I, N);
}

vec3 vec3::refract(const vec3& uv, const vec3& n, double etai_over_etat) {
    auto cos_theta = std::fmin(dot(-uv, n), 1.0);
    vec3 r_out_perp = etai_over_etat * (uv + cos_theta * n);
    vec3 r_out_parallel = -std::sqrt(std::fabs(1.0 - r_out_perp.length_square())) * n;
    return r_out_perp + r_out_parallel;
}

double vec3::random_double() {
    return std::rand() / (RAND_MAX + 1.0);
}

double vec3::random_double(double min, double max) {
    return min + (max - min) * random_double();
}

vec3 vec3::random() {
    return vec3(random_double(), random_double(), random_double());
}

vec3 vec3::random(double min, double max) {
    return vec3(random_double(min, max), random_double(min, max), random_double(min, max));
}

vec3 vec3::random_in_unit_disk() {
    while (true) {
        auto p = vec3(random_double(-1, 1), random_double(-1, 1), 0);
        if (p.length_square() < 1)
            return p;
    }
}

vec3 vec3::random_unit_vector() {
    while (true) {
        auto p = vec3::random(-1, 1);
        auto lensq = p.length_square();
        if (1e-160 < lensq && lensq <= 1.0)
            return p / sqrt(lensq);
    }
}

vec3 vec3::random_in_unit_sphere() {
    while (true) {
        auto p = vec3::random(-1, 1);
        if (p.length_square() >= 1) continue;
        return p;
    }
}

vec3 vec3::random_in_hemisphere(const vec3& normal) {
    vec3 in_unit_sphere = random_in_unit_sphere();
    if (dot(in_unit_sphere, normal) > 0.0) // 同一半球
        return in_unit_sphere;
    else
        return -in_unit_sphere;
}

vec3 vec3::min(const vec3& a, const vec3& b) {
    return vec3(
        std::min(a.x, b.x),
        std::min(a.y, b.y),
        std::min(a.z, b.z)
    );
}

vec3 vec3::max(const vec3& a, const vec3& b) {
    return vec3(
        std::max(a.x, b.x),
        std::max(a.y, b.y),
        std::max(a.z, b.z)
    );
}

bool vec3::is_valid() const {
    return std::isfinite(x) && std::isfinite(y) && std::isfinite(z);
}