// material.h
#ifndef MATERIAL_H
#define MATERIAL_H

#include "vec3.h"
#include "ray.h"
#include "hit_record.h"
#include "interval.h"
#include "photon.h"
#include "photonmap.h"

// �w�q Material ���O�]�����^
class Material {
public:
    virtual ~Material() = default;

    // �������g��ơA�Ҧ��l�����n��{
    virtual bool scatter(
        const Ray& r_in, const HitRecord& rec, vec3& attenuation, Ray& scattered
    ) const = 0;

    virtual vec3 emitted(const Ray& r_in, const HitRecord& rec) const {
        return vec3(0, 0, 0);
    }

    virtual bool photonScatter(
        Photon& photon, const HitRecord& rec, PhotonMap& photonMap
    ) const = 0;

    virtual vec3 directLighting(
        const vec3& normal, const vec3& lightDir, const vec3& lightColor, const vec3& viewDir
    ) const = 0;
};

#endif

