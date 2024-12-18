#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <limits>
#include <algorithm>
#include <random>
#include <unordered_map>
#include <ctime>

#define _USE_MATH_DEFINES
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "vec3.h"
#include "ray.h"
#include "sphere.h"
#include "triangle.h"
#include "material.h"
#include "lambertian.h"
#include "metal.h"
#include "dielectric.h"
#include "photon.h"
#include "photonmap.h"
#include "hit_record.h"
#include "interval.h"
#include "hit_point.h"
#include "hash_grid.h"

const double EPSILON = 1e-4;
const int maxDepth = 5;
const double alpha = 0.7;          // 半徑縮小因子
double maxRadius = 0;

Material* currentMaterial = nullptr;

vec3 lightPosition;
std::vector<Sphere> spheres;
std::vector<Triangle> triangles;

bool findClosestIntersection(const Ray& ray, HitRecord& closest_rec, double t_min, double t_max) {
    HitRecord temp_rec;
    bool hit_anything = false;
    double closest_so_far = t_max;

    // 檢查所有球體
    for (const auto& sphere : spheres) {
        if (sphere.rayIntersect(ray, t_min, closest_so_far, temp_rec)) {
            hit_anything = true;
            closest_so_far = temp_rec.t;
            closest_rec = temp_rec;
        }
    }

    // 檢查所有三角形
    for (const auto& triangle : triangles) {
        if (triangle.rayIntersect(ray, t_min, closest_so_far, temp_rec)) {
            hit_anything = true;
            closest_so_far = temp_rec.t;
            closest_rec = temp_rec;
        }
    }

    return hit_anything;
}

void initializeHitPoints(std::vector<HitPoint>& hitPoints) {
    for (auto& hp : hitPoints) {
        hp.radius2 = initialRadius * initialRadius;
        hp.flux = vec3(0, 0, 0);
        hp.photonCount = 0;
    }
}

void emitPhotonsAndUpdateHitPoints(const vec3& lightPosition, const vec3& lightPower,
    int photonsPerIteration, std::vector<HitPoint>& hitPoints,
    HashGrid& hashGrid) {

    for (int i = 0; i < photonsPerIteration; ++i) {
        // 生成隨機方向
        vec3 dir = vec3::random_in_unit_sphere().normalize();

        // 初始光子
        Photon photon(lightPosition, dir, vec3(20.0, 20.0, 20.0));

        // 追蹤光子
        Ray photonRay(photon.position, photon.direction);

        int depth = 0;
        const int maxPhotonDepth = 5;

        while (depth < maxPhotonDepth) {
            HitRecord rec;

            if (findClosestIntersection(photonRay, rec, EPSILON, std::numeric_limits<double>::max())) {

                if (rec.material_ptr == nullptr) {
                    std::cerr << "Warning: Material pointer is NULL!" << std::endl;
                    return;
                }

                Lambertian* lambertMat = dynamic_cast<Lambertian*>(rec.material_ptr);
                if (lambertMat != nullptr) {
                    // 查找附近的命中點
                    std::vector<size_t> nearbyHitPoints;
                    double queryRadius2 = initialRadius * initialRadius;
                    maxRadius = queryRadius2;
                    hashGrid.query(rec.point, hitPoints[0].radius2, nearbyHitPoints);

                    for (size_t idx : nearbyHitPoints) {
                        HitPoint& hp = hitPoints[idx];
                        double dist2 = (hp.position - rec.point).length_square();
                        if (dist2 <= hp.radius2) {
                            // 累積光通量
                            double M = hp.photonCount + alpha;
                            //double newRadius2 = hp.radius2 * (M / (M + 1.0));
                            //vec3 correctedNormal = vec3::dot(photon.direction, hp.normal) < 0 ? hp.normal : -hp.normal;
                            double cosTerm = vec3::dot(photon.direction, hp.normal) < 0 ? -vec3::dot(photon.direction, hp.normal) : vec3::dot(photon.direction, hp.normal);
                            vec3 phi = hp.flux + hp.throughput * photon.power * cosTerm;

                            hp.photonCount += 1;
                            //hp.radius2 = newRadius2;
                            hp.flux = phi;
                            //hp.flux.print();
                        }
                    }
                    break;
                }

                // 散射光子
                vec3 attenuation;
                if (rec.material_ptr->photonScatter(photon, rec, attenuation)) {
                    photon.power = photon.power * attenuation;
                    //photon.power.print();
                    photonRay = Ray(rec.point + photon.direction * EPSILON, photon.direction);
                    depth++;
                }
                else {
                    break; // 光子被吸收
                }
            }
            else {
                break; // 光子射向無窮遠
            }
        }
    }
}

void traceEyeRay(const Ray& ray, int depth, vec3 throughput, int pixelIndex, std::vector<HitPoint>& hitPoints) {
    if (depth > maxDepth) return;

    HitRecord rec;
    if (findClosestIntersection(ray, rec, 0.001, std::numeric_limits<double>::max())) {

        if (rec.material_ptr == nullptr) {
            std::cerr << "Warning: Material pointer is NULL!" << std::endl;
            return;
        }

        Lambertian* lambertMat = dynamic_cast<Lambertian*>(rec.material_ptr);
        if (lambertMat != nullptr) {
            HitPoint hp;
            hp.position = rec.point;
            hp.normal = rec.normal;
            hp.throughput = throughput * lambertMat->getAlbedo();
            hp.radius2 = initialRadius * initialRadius;
            hp.flux = vec3(0.0, 0.0, 0.0);
            hp.photonCount = 0;
            hp.pixelIndex = pixelIndex;
            hitPoints.push_back(hp);

            return;
        }
        

        // 間接照明
        vec3 attenuation;
        Ray scattered;
        if (rec.material_ptr->scatter(ray, rec, attenuation, scattered)) {
            traceEyeRay(scattered, depth + 1, throughput * attenuation, pixelIndex, hitPoints);
        }
    }
}

Ray computeCameraRay(int x, int y, int width, int height, double fov, const vec3& eyePosition, const vec3& viewDirection, const vec3& upVector) {

    vec3 w = -viewDirection;
    vec3 u = upVector.cross(w).normalize();
    vec3 v = w.cross(u);

    double aspectRatio = static_cast<double>(width) / height;
    double angle = tan(fov * 0.5 * M_PI / 180.0);

    // 將像素座標轉換為NDC
    double ndcX = (x + 0.5) / width;
    double ndcY = (y + 0.5) / height;

    // 將NDC轉換為螢幕空間座標
    double screenX = (2 * ndcX - 1) * aspectRatio * angle;
    double screenY = (1 - 2 * ndcY) * angle;

    vec3 rayDirection = (u * screenX + v * screenY - w).normalize();

    return Ray(eyePosition, rayDirection);
}

// 追蹤視線光線並存儲命中點
void traceEyeRays(const vec3& eyePosition, const vec3& viewDirection, const vec3& upVector,
    double fov, int width, int height, std::vector<HitPoint>& hitPoints,
    std::mt19937& rng, std::uniform_real_distribution<double>& dist) {

    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            // 計算像素對應的視線光線
            Ray ray = computeCameraRay(x, y, width, height, fov, eyePosition, viewDirection, upVector);

            // 追蹤視線光線
            int pixelIndex = y * width + x;
            traceEyeRay(ray, 0, vec3(1.0, 1.0, 1.0), pixelIndex, hitPoints);
        }
    }
}

// 生成最終影像
void generateFinalImage(std::vector<HitPoint>& hitPoints, std::vector<vec3>& image,
    int photonsPerIteration, int numIterations) {
    for (const auto& hp : hitPoints) {
        if (hp.photonCount > 0) {
            // 計算輻射估計
            vec3 radiance = hp.flux / (M_PI * hp.radius2 * photonsPerIteration * numIterations);

            //radiance.print();
            /*
            // 將radiance限制在 [0,1] 範圍內
            radiance.x = std::max(0.0, std::min(1.0, radiance.x));
            radiance.y = std::max(0.0, std::min(1.0, radiance.y));
            radiance.z = std::max(0.0, std::min(1.0, radiance.z));

            // gamma校正
            radiance.x = pow(radiance.x, 1.0 / 2.2);
            radiance.y = pow(radiance.y, 1.0 / 2.2);
            radiance.z = pow(radiance.z, 1.0 / 2.2);
            */

            image[hp.pixelIndex] = image[hp.pixelIndex] + radiance;
            //image[hp.pixelIndex].print();
        }
    }
}

int main(int argc, char* argv[]) {
    // Check input file
    if (argc != 2) {
        std::cerr << "Usage: raytracinghw2 hw2_input.txt" << std::endl;
        return 1;
    }

    std::ifstream infile(argv[1]);
    if (!infile) {
        std::cerr << "Cannot open input file." << std::endl;
        return 1;
    }

    vec3 eyePosition;
    vec3 viewDirection;
    vec3 upVector;
    double fovAngle = 60.0; // Default field of view
    int width = 0, height = 0;

    const int numSamples = 16;

    std::vector<Material*> materials; // 用於管理材質指標
    Material* currentMaterial = nullptr; // 當前材質指標

    std::string line;
    while (std::getline(infile, line)) {
        if (line.empty())
            continue;

        std::istringstream iss(line);
        char type;
        iss >> type;

        if (type == 'E') {
            double x, y, z;
            iss >> x >> y >> z;
            eyePosition = vec3(x, y, z);
        }
        else if (type == 'V') {
            double Dx, Dy, Dz, Ux, Uy, Uz;
            iss >> Dx >> Dy >> Dz >> Ux >> Uy >> Uz;
            viewDirection = vec3(Dx, Dy, Dz).normalize();
            upVector = vec3(Ux, Uy, Uz).normalize();
        }
        else if (type == 'F') {
            iss >> fovAngle; // angle in degrees
        }
        else if (type == 'R') {
            iss >> width >> height;
        }
        else if (type == 'S') {
            double Ox, Oy, Oz, r;
            iss >> Ox >> Oy >> Oz >> r;
            spheres.push_back(Sphere(vec3(Ox, Oy, Oz), r, currentMaterial));
        }
        else if (type == 'T') {
            double x1, y1, z1;
            double x2, y2, z2;
            double x3, y3, z3;
            iss >> x1 >> y1 >> z1
                >> x2 >> y2 >> z2
                >> x3 >> y3 >> z3;
            triangles.push_back(Triangle(vec3(x1, y1, z1),
                vec3(x2, y2, z2),
                vec3(x3, y3, z3),
                currentMaterial));
        }
        else if (type == 'L') {
            double x, y, z;
            iss >> x >> y >> z;
            lightPosition = vec3(x, y, z);
        }
        else if (type == 'M') {
            std::string material_type;
            iss >> material_type;

            if (material_type == "Lambertian") {
                double r, g, b;
                iss >> r >> g >> b;
                currentMaterial = new Lambertian(vec3(r, g, b));
            }
            else if (material_type == "Metal") {
                double r, g, b, fuzz;
                iss >> r >> g >> b >> fuzz;
                currentMaterial = new Metal(vec3(r, g, b), fuzz);
            }
            else if (material_type == "Dielectric") {
                double ref_idx;
                iss >> ref_idx;
                currentMaterial = new Dielectric(ref_idx);
            }

            materials.push_back(currentMaterial); // 將材質指標存入向量
        }
    }

    // 定義光源能量
    vec3 lightPower(30.0, 30.0, 30.0);

    // Check resolution
    if (width == 0 || height == 0) {
        std::cerr << "Invalid resolution." << std::endl;
        return 1;
    }

    // 初始化命中點和影像緩衝區
    std::vector<HitPoint> hitPoints;
    std::vector<vec3> image(width* height, vec3(0, 0, 0));

    // 初始化命中點
    initializeHitPoints(hitPoints);

    // 隨機數生成器
    std::mt19937 rng(static_cast<unsigned int>(time(NULL)));
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    // 追蹤視線光線並存儲命中點
    traceEyeRays(eyePosition, viewDirection, upVector, fovAngle, width, height, hitPoints, rng, dist);

    // 建立hash grid
    HashGrid hashGrid(initialRadius);

    

    // 定義迭代參數
    const int numIterations = 100;              // 迭代次數
    const int photonsPerIteration = 100000;     // 每次迭代發射的光子數
    

    // 迭代 PPM
    for (int iter = 0; iter < numIterations; ++iter) {
        std::cout << "Iteration: " << iter + 1 << "/" << numIterations << std::endl;

        hashGrid.clear();

        // 插入命中點到hash grid
        for (size_t i = 0; i < hitPoints.size(); ++i) {
            if (hitPoints[i].pixelIndex != -1) {
                hashGrid.insert(hitPoints[i], i);
            }
        }

        /*
        int t = 0;
        for (auto& hp : hitPoints) {
            std::cout << "radius_1: " << hp.radius2 << std::endl;
            std::cout << "Photon count_1: " << hp.photonCount << std::endl;
            t++;
            if (t >= 10) {
                break;
            }  
        }
        */
        
        // 發射光子並更新命中點的輻射估計
        emitPhotonsAndUpdateHitPoints(lightPosition, lightPower, photonsPerIteration, hitPoints, hashGrid);

        /*
        t = 0;
        for (auto& hp : hitPoints) {
            std::cout << "radius_2: " << hp.radius2 << std::endl;
            std::cout << "Photon count_2: " << hp.photonCount << std::endl;
            hp.flux.print();
            t++;
            if (t >= 10) {
                break;
            }
        }
        */

        // 縮小命中點的搜索半徑
        for (auto& hp : hitPoints) {
            if (hp.photonCount > 0) {
                hp.radius2 *= (static_cast<double>(hp.photonCount) + alpha) / (static_cast<double>(hp.photonCount) + 1.0);
                hp.flux *= (static_cast<double>(hp.photonCount) + alpha) / (static_cast<double>(hp.photonCount) + 1.0);

                if (hp.radius2 > maxRadius) {
                    maxRadius = hp.radius2;
                }
            }
        }
    }

    // 生成最終影像
    generateFinalImage(hitPoints, image, photonsPerIteration, numIterations);

    // 保存影像為 PPM 檔案
    std::string filename = "output.ppm";
    std::ofstream ppmFile(filename, std::ios::binary);
    ppmFile << "P6\n" << width << " " << height << "\n255\n";
    for (const auto& color : image) {
        unsigned char r = static_cast<unsigned char>(std::min(color.x * 255.0, 255.0));
        unsigned char g = static_cast<unsigned char>(std::min(color.y * 255.0, 255.0));
        unsigned char b = static_cast<unsigned char>(std::min(color.z * 255.0, 255.0));
        ppmFile.write(reinterpret_cast<char*>(&r), 1);
        ppmFile.write(reinterpret_cast<char*>(&g), 1);
        ppmFile.write(reinterpret_cast<char*>(&b), 1);
    }

    ppmFile.close();    

    std::cout << "Image has been saved to " << filename << std::endl;

    // 釋放材質指標
    for (auto material_ptr : materials) {
        delete material_ptr;
    }

    return 0;
}

