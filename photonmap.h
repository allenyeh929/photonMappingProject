#ifndef PHOTONMAP_H
#define PHOTONMAP_H

#include <vector>
#include <queue>
#include <functional>
#include "photon.h"
#include "vec3.h"

struct KDTreeNode {
    Photon photon;
    KDTreeNode* left;
    KDTreeNode* right;
    int plane; // ���ζb�A0: x, 1: y, 2: z

    KDTreeNode(const Photon& p)
        : photon(p), left(nullptr), right(nullptr), plane(0) {}
};

struct PhotonDist {
    const Photon* photon;
    double dist2;

    PhotonDist(const Photon* p, double d) : photon(p), dist2(d) {}

    // �w�q����B��l�A�إ̤߳j��
    bool operator<(const PhotonDist& other) const {
        return dist2 < other.dist2;
    }
};

class PhotonMap {
public:
    KDTreeNode* root;

    PhotonMap();
    ~PhotonMap();

    void store(const Photon& photon); // �x�s���l
    void balance();                   // ����kd-tree
    void locatePhotons(const vec3& position, int maxPhotons,
        std::vector<const Photon*>& foundPhotons) const; // �d�ߥ��l

private:
    std::vector<Photon> photons;

    KDTreeNode* buildKDTree(std::vector<Photon>& photons, int start, int end);
    void deleteKDTree(KDTreeNode* node);
    void locatePhotonsRecursive(const KDTreeNode* node,
        const vec3& position, int maxPhotons,
        std::priority_queue<PhotonDist>& photonHeap) const;
};

#endif

