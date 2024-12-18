#ifndef HIT_POINT_H
#define HIT_POINT_H

#include "vec3.h"

const double initialRadius = 0.2;

struct HitPoint {
    vec3 position;      // 命中點的位置
    vec3 normal;        // 命中點的法線
    vec3 flux;          // 累積的光通量
    vec3 throughput;    // 從視線光線到此點的輻射累積
    double radius2;     // 當前的搜索半徑平方
    int photonCount;    // 累積的光子數
    int pixelIndex;     // 影像緩衝區中的像素索引

    HitPoint()
        : position(0, 0, 0), normal(0, 0, 0), flux(0, 0, 0),
        throughput(1, 1, 1), radius2(initialRadius * initialRadius),
        photonCount(0), pixelIndex(-1) {}
};

#endif

