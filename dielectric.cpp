#include "dielectric.h"

bool Dielectric::scatter(
    const Ray& r_in, const HitRecord& rec, vec3& attenuation, Ray& scattered
) const {
    attenuation = vec3(1.0, 1.0, 1.0); // �������褣�l������q
    double etai_over_etat = rec.front_face ? (1.0 / ref_idx) : ref_idx;

    vec3 unit_direction = r_in.direction.normalize();
    double cos_theta = fmin(vec3::dot(-unit_direction, rec.normal), 1.0);
    double sin_theta = sqrt(1.0 - cos_theta * cos_theta);

    bool cannot_refract = etai_over_etat * sin_theta > 1.0;
    vec3 direction;

    /*
    if (cannot_refract || reflectance(cos_theta, etai_over_etat) > random_double(0.0, 1.0)) {
        // �Ϯg
        //direction = vec3::reflect(unit_direction, rec.normal);
    }
    else {
        // ��g
        direction = vec3::refract(unit_direction, rec.normal, etai_over_etat);
    }
    */

    direction = vec3::refract(unit_direction, rec.normal, etai_over_etat);

    direction = direction.normalize();
    scattered = Ray(rec.point, direction);
    return true;
}

bool Dielectric::photonScatter(
    Photon& photon, const HitRecord& rec, vec3& attenuation
) const {
    vec3 unit_direction = photon.direction;
    double cos_theta = fmin(vec3::dot(-unit_direction, rec.normal), 1.0);
    double sin_theta = sqrt(1.0 - cos_theta * cos_theta);

    double etai_over_etat = rec.front_face ? (1.0 / ref_idx) : ref_idx;
    bool cannot_refract = etai_over_etat * sin_theta > 1.0;

    vec3 direction;

    /*
    if (cannot_refract || reflectance(cos_theta, etai_over_etat) > random_double(0.0, 1.0)) {
        // �Ϯg
        //direction = vec3::reflect(unit_direction, rec.normal);
    }
    else {
        // ��g
        direction = vec3::refract(unit_direction, rec.normal, etai_over_etat);
    }
    */

    direction = vec3::refract(unit_direction, rec.normal, etai_over_etat);

    direction = direction.normalize();
    photon.direction = direction;
    attenuation = vec3(1.0, 1.0, 1.0);

    return true; // ���l�~��Ǽ�
}

vec3 Dielectric::directLighting(
    const vec3& normal, const vec3& lightDir, const vec3& lightColor, const vec3& viewDir
) const {
    return vec3(0, 0, 0); // ���Ҽ{��������
}