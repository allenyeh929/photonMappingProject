#include "metal.h"

bool Metal::scatter(
    const Ray& r_in, const HitRecord& rec, vec3& attenuation, Ray& scattered
) const {
    vec3 reflected = vec3::reflect(r_in.direction.normalize(), rec.normal);
    reflected = reflected + fuzz * vec3::random_in_unit_sphere();
    scattered = Ray(rec.point, reflected);
    attenuation = albedo;
    return (vec3::dot(scattered.direction, rec.normal) > 0);
}

bool Metal::photonScatter(
    Photon& photon, const HitRecord& rec, vec3& attenuation

) const {
    // Russian Roulette 決定光子是否被吸收
    double survivalProbability = albedo.length() / sqrt(3);
    if (random_double(0.0, 1.0) > survivalProbability) {
        return false; // 光子被吸收
    }

    // 更新光子的能量
    photon.power = photon.power * albedo / survivalProbability;

    // 鏡面反射新的方向
    photon.direction = vec3::reflect(photon.direction, rec.normal);
    photon.direction = photon.direction + fuzz * vec3::random_in_unit_sphere();
    photon.direction = photon.direction.normalize();

    return (vec3::dot(photon.direction, rec.normal) > 0); // 光子是否繼續傳播
}

vec3 Metal::directLighting(
    const vec3& normal, const vec3& lightDir, const vec3& lightColor, const vec3& viewDir
) const {
    vec3 reflected = vec3::reflect(-lightDir, normal).normalize();
    double specularFactor = pow(std::max(0.0, vec3::dot(reflected, viewDir)), 50); // 50 為高光指數
    return albedo * lightColor * specularFactor;
}