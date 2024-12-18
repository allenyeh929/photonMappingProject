#ifndef HIT_POINT_H
#define HIT_POINT_H

#include "vec3.h"

const double initialRadius = 0.2;

struct HitPoint {
    vec3 position;      // �R���I����m
    vec3 normal;        // �R���I���k�u
    vec3 flux;          // �ֿn�����q�q
    vec3 throughput;    // �q���u���u�즹�I����g�ֿn
    double radius2;     // ��e���j���b�|����
    int photonCount;    // �ֿn�����l��
    int pixelIndex;     // �v���w�İϤ�����������

    HitPoint()
        : position(0, 0, 0), normal(0, 0, 0), flux(0, 0, 0),
        throughput(1, 1, 1), radius2(initialRadius * initialRadius),
        photonCount(0), pixelIndex(-1) {}
};

#endif

