#ifndef PHOTON_H
#define PHOTON_H

#include "vec3.h"

class Photon {
public:
	vec3 position; // 光子的位置
	vec3 direction; // 光子的入射方向
	vec3 power; // 光子的能量(顏色)
	int plane; // 分割軸（0: x, 1: y, 2: z）

	Photon();
	Photon(const vec3& position_, const vec3& direction_, const vec3& power_);
};

#endif

