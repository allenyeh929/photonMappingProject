#include "photonmap.h"
#include <algorithm>
#include <iostream>
#include <cmath>

PhotonMap::PhotonMap() : root(nullptr) {}

PhotonMap::~PhotonMap() {
    deleteKDTree(root);
}

void PhotonMap::deleteKDTree(KDTreeNode* node) {
    if (node) {
        deleteKDTree(node->left);
        deleteKDTree(node->right);
        delete node;
    }
}

void PhotonMap::store(const Photon& photon) {
    if (!photon.position.is_valid() || !photon.direction.is_valid()) {
        std::cerr << "Invalid photon detected, skipping store." << std::endl;
        return;
    }
    photons.push_back(photon);
}

void PhotonMap::balance() {
    if (photons.empty()) return;
    root = buildKDTree(photons, 0, photons.size());
    photons.clear(); // 清除暫存的光子數據
}

KDTreeNode* PhotonMap::buildKDTree(std::vector<Photon>& photons,
    int start, int end) {
    if (start >= end) return nullptr;

    // 計算 bounding box，決定分割軸
    vec3 minBound = photons[start].position;
    vec3 maxBound = photons[start].position;
    for (int i = start + 1; i < end; ++i) {
        minBound = vec3::min(minBound, photons[i].position);
        maxBound = vec3::max(maxBound, photons[i].position);
    }

    int axis = 0;
    vec3 extent = maxBound - minBound;
    if (extent.y > extent.x && extent.y > extent.z)
        axis = 1;
    else if (extent.z > extent.x)
        axis = 2;

    // 按照分割軸對光子排序，找到中位數
    int mid = (start + end) / 2;
    std::nth_element(photons.begin() + start, photons.begin() + mid,
        photons.begin() + end,
        [axis](const Photon& a, const Photon& b) {
            return a.position[axis] < b.position[axis];
        });

    // 創建節點
    KDTreeNode* node = new KDTreeNode(photons[mid]);
    node->plane = axis;

    // 遞迴構建子樹
    node->left = buildKDTree(photons, start, mid);
    node->right = buildKDTree(photons, mid + 1, end);

    return node;
}

void PhotonMap::locatePhotons(const vec3& position, int maxPhotons,
    std::vector<const Photon*>& foundPhotons) const {
    
    std::priority_queue<PhotonDist> photonHeap;

    // 呼叫遞迴函數
    locatePhotonsRecursive(root, position, maxPhotons, photonHeap);

    // 將結果從堆中取出
    foundPhotons.clear();
    while (!photonHeap.empty()) {
        foundPhotons.push_back(photonHeap.top().photon);
        photonHeap.pop();
    }
}

void PhotonMap::locatePhotonsRecursive(const KDTreeNode* node,
    const vec3& position, int maxPhotons,
    std::priority_queue<PhotonDist>& photonHeap) const {
    if (!node) return;

    // 計算當前節點與查詢點的距離平方
    double dist2 = (node->photon.position - position).length_square();

    // 更新priority_queue
    if (photonHeap.size() < maxPhotons) {
        photonHeap.emplace(&node->photon, dist2);
    }
    else if (dist2 < photonHeap.top().dist2) {
        photonHeap.pop();
        photonHeap.emplace(&node->photon, dist2);
    }

    // 遞迴遍歷子樹
    int axis = node->plane;
    double diff = position[axis] - node->photon.position[axis];

    KDTreeNode* nearChild = diff < 0 ? node->left : node->right;
    KDTreeNode* farChild = diff < 0 ? node->right : node->left;

    // 先遍歷與查詢點較近的一側
    locatePhotonsRecursive(nearChild, position, maxPhotons, photonHeap);

    // 判斷是否需要遍歷另一側
    if (photonHeap.size() < maxPhotons || diff * diff < photonHeap.top().dist2) {
        locatePhotonsRecursive(farChild, position, maxPhotons, photonHeap);
    }
}