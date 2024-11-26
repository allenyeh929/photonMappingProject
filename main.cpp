#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <limits>
#include <algorithm>

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

const double EPSILON = 1e-4;
const int maxDepth = 5;

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


void emitPhotons(PhotonMap& photonMap, const vec3& lightPosition, const vec3& lightPower, int numPhotons) {
    for (int i = 0; i < numPhotons; ++i) {
        // 生成隨機方向
        vec3 dir = vec3::random_in_unit_sphere().normalize();

        // 初始光子
        Photon photon(lightPosition, dir, lightPower / numPhotons);

        // 追蹤光子
        Ray photonRay(photon.position, photon.direction);

        int depth = 0;
        const int maxPhotonDepth = 5;

        while (depth < maxPhotonDepth) {
            HitRecord rec;
            if (findClosestIntersection(photonRay, rec, EPSILON, std::numeric_limits<double>::max())) {
                // 獲取交點處的材質
                Material* material = rec.material_ptr;

                if (!material->photonScatter(photon, rec, photonMap)) {
                    break; // 光子被吸收
                }

                // 更新光子的射線
                photonRay = Ray(photon.position + photon.direction * EPSILON, photon.direction);

                depth++;
            }
            else {
                break; // 光子射向無窮遠
            }
        }
    }

    photonMap.balance(); // 平衡kd-tree
}

bool isInShadow(const vec3& point, const vec3& lightDir) {
    Ray shadowRay(point + lightDir * EPSILON, lightDir);
    double lightDistance = (lightPosition - point).length();

    HitRecord rec;

    // 檢查與場景中所有物體的交點
    if (findClosestIntersection(shadowRay, rec, EPSILON, lightDistance)) {
        return true; // 有遮擋
    }

    return false; // 無遮擋
}

vec3 computeLighting(const vec3& point, const vec3& normal, const vec3& viewDir, const Material* material, const PhotonMap& photonMap) {
    // Ambient component
    vec3 ambient = vec3(0, 0, 0);

    // 直接光照
    vec3 lightDir = (lightPosition - point).normalize();
    vec3 diffuse(0, 0, 0);
    vec3 specular(0, 0, 0);

    vec3 lightColor(1, 1, 1);
    // 檢查陰影（判斷光線是否被遮擋）
    if (!isInShadow(point, lightDir)) {
        // 調用材質的 directLighting 函數
        diffuse = material->directLighting(normal, lightDir, lightColor, viewDir);
    }

    // 間接光照（使用光子映射）
    vec3 indirect(0, 0, 0);
    const int maxPhotons = 50;
    std::vector<const Photon*> photons;
    photonMap.locatePhotons(point, maxPhotons, photons);

    if (!photons.empty()) {
        //std::cout << "ok" << std::endl;
        double maxDist2 = 0.0;
        for (const Photon* photon : photons) {
            double dist2 = (photon->position - point).length_square();
            if (dist2 > maxDist2) {
                maxDist2 = dist2;
            }
        }
        
        double area = M_PI * maxDist2;

        vec3 flux(0, 0, 0);
        for (const Photon* photon : photons) {
            flux = flux + photon->power;
        }

        indirect = flux / area;
    }

    // 合併光照
    vec3 color = ambient + diffuse + specular + indirect;

    color = color / (color + vec3(1.0, 1.0, 1.0));

    // 伽馬校正
    color.x = pow(color.x, 1.0 / 2.2);
    color.y = pow(color.y, 1.0 / 2.2);
    color.z = pow(color.z, 1.0 / 2.2);

    /*
    // 限制顏色值在 [0,1] 範圍內，防止過曝
    color.x = std::min(1.0, std::max(0.0, color.x));
    color.y = std::min(1.0, std::max(0.0, color.y));
    color.z = std::min(1.0, std::max(0.0, color.z));
    */

    // 返回計算後的顏色
    return color;
}

vec3 trace(const Ray& ray, int depth, const PhotonMap& photonMap) {
    if (depth > maxDepth) {
        return vec3(0, 0, 0); // Black
    }

    HitRecord rec;
    if (findClosestIntersection(ray, rec, 0.001, std::numeric_limits<double>::max())) {
        Ray scattered;
        vec3 attenuation;
        vec3 emitted = rec.material_ptr->emitted(ray, rec);

        if (rec.material_ptr->scatter(ray, rec, attenuation, scattered)) {
            vec3 indirectLighting = trace(scattered, depth + 1, photonMap);
            vec3 directLighting = computeLighting(rec.point, rec.normal, -ray.direction, rec.material_ptr, photonMap);

            return emitted + attenuation * 0.8 * (indirectLighting + directLighting);
        }
        else {
            return emitted;
        }
    }
    else {
        // 背景色
        return vec3(0, 0, 0);
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

    const int numSamples = 32;

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

    // 定義光子映射
    PhotonMap photonMap;

    // 定義光源能量
    vec3 lightPower(0.1, 0.1, 0.1);

    // 發射光子
    emitPhotons(photonMap, lightPosition, lightPower, 1000000);

    // Check resolution
    if (width == 0 || height == 0) {
        std::cerr << "Invalid resolution." << std::endl;
        return 1;
    }

    // Prepare image buffer
    std::vector<unsigned char> image(width * height * 3, 0); // Initialize to black

    // Compute camera basis vectors
    double aspectRatio = static_cast<double>(width) / height;
    vec3 right = upVector.cross(viewDirection).normalize();
    vec3 up = viewDirection.cross(right).normalize();

    // Viewing distance (arbitrary)
    double d = 1.0;

    // Compute image plane dimensions
    double fovRadians = fovAngle * M_PI / 180.0;
    double halfWidth = tan(fovRadians / 2.0);
    double halfHeight = halfWidth / aspectRatio;

    // For each pixel
    for (int j = 0; j < height; ++j) {
        for (int i = 0; i < width; ++i) {
            vec3 color = vec3(0, 0, 0);

            for (int s = 0; s < numSamples; ++s) {
                // Compute u and v
                double u = ((i + 0.5) / width - 0.5) * 2 * halfWidth;
                double v = (0.5 - (j + 0.5) / height) * 2 * halfHeight;

                // Compute ray direction
                vec3 rayDirection = (viewDirection * d + right * u + up * v).normalize();

                Ray ray(eyePosition, rayDirection);

                // Trace the ray
                color = color + trace(ray, 0, photonMap);
            }

            // Average the color samples
            color = color * (1.0 / numSamples);

            // Convert color to [0,255] and write to image
            int index = (j * width + i) * 3;
            if (index + 2 >= image.size()) {
                std::cerr << "Pixel index out of range: " << index << std::endl;
                return 1;
            }

            image[index] = static_cast<unsigned char>(std::min(color.x * 255.0, 255.0));
            image[index + 1] = static_cast<unsigned char>(std::min(color.y * 255.0, 255.0));
            image[index + 2] = static_cast<unsigned char>(std::min(color.z * 255.0, 255.0));
        }
    }

    // Save the image as a PPM file
    std::string filename = "output.ppm";
    std::ofstream ppmFile(filename, std::ios::binary);
    ppmFile << "P6\n" << width << " " << height << "\n255\n";
    ppmFile.write(reinterpret_cast<char*>(image.data()), image.size());
    ppmFile.close();    

    std::cout << "Image has been saved to " << filename << std::endl;

    // 釋放材質指標
    for (auto material_ptr : materials) {
        delete material_ptr;
    }

    return 0;
}
