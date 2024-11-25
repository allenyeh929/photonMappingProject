#include "dielectric.h"

bool Dielectric::scatter(
    const Ray& r_in, const HitRecord& rec, vec3& attenuation, Ray& scattered
) const {
    attenuation = vec3(1.0, 1.0, 1.0); // 玻璃材質不吸收光能量
    double etai_over_etat = rec.front_face ? (1.0 / ref_idx) : ref_idx;

    vec3 unit_direction = r_in.direction.normalize();
    double cos_theta = fmin(vec3::dot(-unit_direction, rec.normal), 1.0);
    double sin_theta = sqrt(1.0 - cos_theta * cos_theta);

    bool cannot_refract = etai_over_etat * sin_theta > 1.0;
    vec3 direction;

    if (cannot_refract || reflectance(cos_theta, etai_over_etat) > random_double(0.0, 1.0)) {
        // 反射
        direction = vec3::reflect(unit_direction, rec.normal);
    }
    else {
        // 折射
        direction = vec3::refract(unit_direction, rec.normal, etai_over_etat);
    }

    scattered = Ray(rec.point, direction);
    return true;
}

bool Dielectric::photonScatter(
    Photon& photon, const HitRecord& rec, PhotonMap& photonMap
) const {
    // 玻璃材質假設不吸收光子
    photon.position = rec.point;
    photonMap.store(photon);

    vec3 unit_direction = photon.direction;
    double cos_theta = fmin(vec3::dot(-unit_direction, rec.normal), 1.0);
    double sin_theta = sqrt(1.0 - cos_theta * cos_theta);

    double etai_over_etat = rec.front_face ? (1.0 / ref_idx) : ref_idx;
    bool cannot_refract = etai_over_etat * sin_theta > 1.0;

    vec3 direction;

    if (cannot_refract || reflectance(cos_theta, etai_over_etat) > random_double(0.0, 1.0)) {
        // 反射
        direction = vec3::reflect(unit_direction, rec.normal);
    }
    else {
        // 折射
        direction = vec3::refract(unit_direction, rec.normal, etai_over_etat);
    }

    photon.direction = direction.normalize();

    return true; // 光子繼續傳播
}

vec3 Dielectric::directLighting(
    const vec3& normal, const vec3& lightDir, const vec3& lightColor, const vec3& viewDir
) const {
    // 對於透明材質，直接光照可忽略或進行更複雜的計算
    return vec3(0, 0, 0); // 返回零向量，表示不考慮直接光照
}