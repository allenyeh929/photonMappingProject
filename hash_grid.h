// hash_grid.h
#ifndef HASH_GRID_H
#define HASH_GRID_H

#include "vec3.h"
#include "hit_point.h"

#include <unordered_map>
#include <vector>
#include <tuple>

// 定義網格座標
struct GridCoord {
    int x, y, z;

    GridCoord() : x(0), y(0), z(0) {}

    GridCoord(int x_, int y_, int z_) : x(x_), y(y_), z(z_) {}

    bool operator==(const GridCoord& other) const {
        return x == other.x && y == other.y && z == other.z;
    }
};

// hash function
struct GridCoordHash {
    std::size_t operator()(const GridCoord& coord) const {
        return std::hash<int>()(coord.x) ^ (std::hash<int>()(coord.y) << 1) ^ (std::hash<int>()(coord.z) << 2);
    }
};

// 定義grid cell
struct GridCell {
    std::vector<size_t> hitPointIndices; // 命中點向量中的索引
};

// 定義hash grid class
class HashGrid {
public:
    HashGrid(double cellSize_) : cellSize(cellSize_) {}

    void clear() {
        grid.clear();
    }

    void insert(const HitPoint& hp, size_t index) {
        GridCoord coord = getGridCoord(hp.position);
        grid[coord].hitPointIndices.push_back(index);
    }

    // 查詢附近的命中點索引
    void query(const vec3& position, double radius2, std::vector<size_t>& results) const {
        int ix = static_cast<int>(floor(position.x / cellSize));
        int iy = static_cast<int>(floor(position.y / cellSize));
        int iz = static_cast<int>(floor(position.z / cellSize));

        int delta = static_cast<int>(ceil(sqrt(radius2) / cellSize));

        for (int dx = -delta; dx <= delta; ++dx) {
            for (int dy = -delta; dy <= delta; ++dy) {
                for (int dz = -delta; dz <= delta; ++dz) {
                    GridCoord coord(ix + dx, iy + dy, iz + dz);
                    auto it = grid.find(coord);
                    if (it != grid.end()) {
                        for (size_t idx : it->second.hitPointIndices) {
                            results.push_back(idx);
                        }
                    }
                }
            }
        }
    }

private:
    double cellSize;
    std::unordered_map<GridCoord, GridCell, GridCoordHash> grid;

    GridCoord getGridCoord(const vec3& position) const {
        return GridCoord{
            static_cast<int>(floor(position.x / cellSize)),
            static_cast<int>(floor(position.y / cellSize)),
            static_cast<int>(floor(position.z / cellSize))
        };
    }
};

#endif

