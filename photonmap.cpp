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

void PhotonMap::locatePhotons(const vec3& position, int maxPhotons,
    std::vector<const Photon*>& foundPhotons) const {
    
    std::priority_queue<PhotonDist> photonHeap;

    // �I�s���j���
    locatePhotonsRecursive(root, position, maxPhotons, photonHeap);

    // �N���G�q�襤���X
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

    // �p���e�`�I�P�d���I���Z������
    double dist2 = (node->photon.position - position).length_square();

    // ��spriority_queue
    if (photonHeap.size() < maxPhotons) {
        photonHeap.emplace(&node->photon, dist2);
    }
    else if (dist2 < photonHeap.top().dist2) {
        photonHeap.pop();
        photonHeap.emplace(&node->photon, dist2);
    }

    // ���j�M���l��
    int axis = node->plane;
    double diff = position[axis] - node->photon.position[axis];

    KDTreeNode* nearChild = diff < 0 ? node->left : node->right;
    KDTreeNode* farChild = diff < 0 ? node->right : node->left;

    // ���M���P�d���I���񪺤@��
    locatePhotonsRecursive(nearChild, position, maxPhotons, photonHeap);

    // �P�_�O�_�ݭn�M���t�@��
    if (photonHeap.size() < maxPhotons || diff * diff < photonHeap.top().dist2) {
        locatePhotonsRecursive(farChild, position, maxPhotons, photonHeap);
    }
}