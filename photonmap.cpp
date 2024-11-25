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
    photons.clear(); // �M���Ȧs�����l�ƾ�
}

KDTreeNode* PhotonMap::buildKDTree(std::vector<Photon>& photons,
    int start, int end) {
    if (start >= end) return nullptr;

    // �p�� bounding box�A�M�w���ζb
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

    // ���Ӥ��ζb����l�ƧǡA��줤���
    int mid = (start + end) / 2;
    std::nth_element(photons.begin() + start, photons.begin() + mid,
        photons.begin() + end,
        [axis](const Photon& a, const Photon& b) {
            return a.position[axis] < b.position[axis];
        });

    // �Ыظ`�I
    KDTreeNode* node = new KDTreeNode(photons[mid]);
    node->plane = axis;

    // ���j�c�ؤl��
    node->left = buildKDTree(photons, start, mid);
    node->right = buildKDTree(photons, mid + 1, end);

    return node;
}

void PhotonMap::locatePhotons(const vec3& position, double maxDist,
    int maxPhotons,
    std::vector<const Photon*>& foundPhotons) const {
    foundPhotons.clear();
    locatePhotonsRecursive(root, position, maxDist * maxDist, maxPhotons,
        foundPhotons);
}

void PhotonMap::locatePhotonsRecursive(const KDTreeNode* node,
    const vec3& position, double maxDist2,
    int maxPhotons,
    std::vector<const Photon*>& foundPhotons) const {
    if (!node) return;

    // �p���e�`�I�P�d���I���Z������
    double dist2 = (node->photon.position - position).length_square();
    if (dist2 < maxDist2) {
        if (foundPhotons.size() < maxPhotons) {
            foundPhotons.push_back(&node->photon);
        }
    }

    // ���j�M���l��
    int axis = node->plane;
    double diff = position[axis] - node->photon.position[axis];

    if (diff < 0) {
        locatePhotonsRecursive(node->left, position, maxDist2,
            maxPhotons, foundPhotons);
        if (diff * diff < maxDist2) {
            locatePhotonsRecursive(node->right, position, maxDist2,
                maxPhotons, foundPhotons);
        }
    }
    else {
        locatePhotonsRecursive(node->right, position, maxDist2,
            maxPhotons, foundPhotons);
        if (diff * diff < maxDist2) {
            locatePhotonsRecursive(node->left, position, maxDist2,
                maxPhotons, foundPhotons);
        }
    }
}

