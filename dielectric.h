#ifndef DIELECTRIC_H
#define DIELECTRIC_H

#include "material.h"
#include "vec3.h"
#include "photon.h"
#include "photonmap.h"

class Dielectric : public Material {
public:
    Dielectric(double refractive_index) : ref_idx(refractive_index) {}

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
    double ref_idx;

    static double reflectance(double cosine, double ref_idx) {
        // Schlick's approximation
        auto r0 = (1 - ref_idx) / (1 + ref_idx);
        r0 = r0 * r0;
        return r0 + (1 - r0) * pow((1 - cosine), 5);
    }
};

#endif

