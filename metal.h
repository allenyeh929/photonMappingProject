#ifndef METAL_H
#define METAL_H

#include "material.h"
#include "vec3.h"
#include "photon.h"
#include "photonmap.h"

class Metal : public Material {
public:
    Metal(const vec3& albedo, double fuzziness) : albedo(albedo), fuzz(fuzziness < 1 ? fuzziness : 1) {}

    virtual bool scatter(
        const Ray& r_in, const HitRecord& rec, vec3& attenuation, Ray& scattered
    ) const override;

    virtual bool photonScatter(
        Photon& photon, const HitRecord& rec, PhotonMap& photonMap
    ) const override;

    virtual vec3 directLighting(
        const vec3& normal, const vec3& lightDir, const vec3& lightColor, const vec3& viewDir
    ) const override;

private:
    vec3 albedo;
    double fuzz;
};

#endif

