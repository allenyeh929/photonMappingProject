#ifndef LAMBERTIAN_H
#define LAMBERTIAN_H

#include "material.h"

class Lambertian : public Material {
public:
    Lambertian(const vec3& albedo) : albedo(albedo) {}

    vec3 getAlbedo() const {
        return albedo;
    }

    virtual bool scatter(
        const Ray& r_in, const HitRecord& rec, vec3& attenuation, Ray& scattered
    ) const override;

    virtual bool photonScatter(
        Photon& photon, const HitRecord& rec, vec3& attenuation
    ) const override;

    virtual vec3 directLighting(
        const vec3& normal, const vec3& lightDir, const vec3& lightColor, const vec3& viewDir
    ) const override;

private:
    vec3 albedo;
};

#endif

