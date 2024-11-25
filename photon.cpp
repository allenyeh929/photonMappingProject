#include "photon.h"

Photon::Photon() : plane(0) {}

Photon::Photon(const vec3& position_, const vec3& direction_, const vec3& power_)
    : position(position_), direction(direction_), power(power_), plane(0) {}