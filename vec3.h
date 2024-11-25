#ifndef VEC3_H
#define VEC3_H

class vec3 {
public:
    double x, y, z;

    vec3();
    vec3(double x_, double y_, double z_);

    double& operator[](int index);
    const double& operator[](int index) const;

    vec3 operator+(const vec3& v) const;
    vec3 operator-(const vec3& v) const;
    vec3 operator*(double s) const;
    vec3 operator/(double s) const;

    vec3 operator-() const;

    vec3 operator*(const vec3& other) const; // 向量與向量的逐元素相乘
    vec3& operator*=(double scalar);            // 向量與純量的乘法賦值

    vec3 operator/(const vec3& other) const;


    static double dot(const vec3& u, const vec3& v);
    vec3 cross(const vec3& v) const;

    double length() const;
    double length_square() const;

    bool near_zero() const;

    vec3 normalize() const;

    static vec3 reflect(const vec3& I, const vec3& N);
    static vec3 refract(const vec3& uv, const vec3& n, double etai_over_etat);

    static double random_double();
    static double random_double(double min, double max);

    static vec3 random();
    static vec3 random(double min, double max);

    static vec3 random_in_unit_disk();
    static vec3 random_unit_vector();

    static vec3 random_in_unit_sphere();
    static vec3 random_in_hemisphere(const vec3& normal);

    static vec3 min(const vec3& a, const vec3& b);
    static vec3 max(const vec3& a, const vec3& b);

    bool is_valid() const;
};

inline vec3 operator*(double t, const vec3& v) {
    return vec3(t * v.x, t * v.y, t * v.z);
}

#endif


