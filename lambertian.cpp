#include "lambertian.h"

bool Lambertian::scatter(
    const Ray& r_in, const HitRecord& rec, vec3& attenuation, Ray& scattered
) const {
    vec3 scatter_direction = rec.normal + vec3::random_unit_vector();

    // Catch degenerate scatter direction
    if (scatter_direction.near_zero())
        scatter_direction = rec.normal;

    scattered = Ray(rec.point, scatter_direction);
    attenuation = albedo;
    return true;
}

bool Lambertian::photonScatter(
    Photon& photon, const HitRecord& rec, PhotonMap& photonMap
) const {
    // �x�s���l
    photon.position = rec.point;
    photonMap.store(photon);

    double survivalProbability = albedo.length() / sqrt(3);
    if (random_double(0.0, 1.0) > survivalProbability) {
        return false; // ���l�Q�l��
    }

    // ��s���l����q
    photon.power = photon.power * albedo / survivalProbability;

    // ���Ϯg�s����V
    photon.direction = vec3::random_in_hemisphere(rec.normal).normalize();

    return true;
}

vec3 Lambertian::directLighting(
    const vec3& normal, const vec3& lightDir, const vec3& lightColor, const vec3& viewDir
) const {
    double NdotL = std::max(0.0, vec3::dot(normal, lightDir));
    return albedo * lightColor * NdotL;
}