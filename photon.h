#ifndef PHOTON_H
#define PHOTON_H

#include "vec3.h"

class Photon {
public:
	vec3 position; // ���l����m
	vec3 direction; // ���l���J�g��V
	vec3 power; // ���l����q(�C��)
	int plane; // ���ζb�]0: x, 1: y, 2: z�^

	Photon();
	Photon(const vec3& position_, const vec3& direction_, const vec3& power_);
};

#endif

