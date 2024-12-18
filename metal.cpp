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
    // Russian Roulette �M�w���l�O�_�Q�l��
    double survivalProbability = albedo.length() / sqrt(3);
    if (random_double(0.0, 1.0) > survivalProbability) {
        return false; // ���l�Q�l��
    }

    // ��s���l����q
    photon.power = photon.power * albedo / survivalProbability;

    // �譱�Ϯg�s����V
    photon.direction = vec3::reflect(photon.direction, rec.normal);
    photon.direction = photon.direction + fuzz * vec3::random_in_unit_sphere();
    photon.direction = photon.direction.normalize();

    return (vec3::dot(photon.direction, rec.normal) > 0); // ���l�O�_�~��Ǽ�
}

vec3 Metal::directLighting(
    const vec3& normal, const vec3& lightDir, const vec3& lightColor, const vec3& viewDir
) const {
    vec3 reflected = vec3::reflect(-lightDir, normal).normalize();
    double specularFactor = pow(std::max(0.0, vec3::dot(reflected, viewDir)), 50); // 50 ����������
    return albedo * lightColor * specularFactor;
}