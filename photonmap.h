#ifndef PHOTONMAP_H
#define PHOTONMAP_H

#include <vector>
#include "photon.h"
#include "vec3.h"

struct KDTreeNode {
    Photon photon;
    KDTreeNode* left;
    KDTreeNode* right;
    int plane; // 分割軸，0: x, 1: y, 2: z

    KDTreeNode(const Photon& p)
        : photon(p), left(nullptr), right(nullptr), plane(0) {}
};

class PhotonMap {
public:
    KDTreeNode* root;

    PhotonMap();
    ~PhotonMap();

    void store(const Photon& photon); // 儲存光子
    void balance();                   // 平衡kd-tree
    void locatePhotons(const vec3& position, double maxDist, int maxPhotons, std::vector<const Photon*>& foundPhotons) const; // 查詢光子

private:
    std::vector<Photon> photons;

    KDTreeNode* buildKDTree(std::vector<Photon>& photons, int start, int end);
    void deleteKDTree(KDTreeNode* node);
    void locatePhotonsRecursive(const KDTreeNode* node, const vec3& position,
        double maxDist2, int maxPhotons,
        std::vector<const Photon*>& foundPhotons) const;
};

#endif // PHOTONMAP_H

