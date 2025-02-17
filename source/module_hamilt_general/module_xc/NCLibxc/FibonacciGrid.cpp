#include <iostream>
#include <vector>
#include <cmath>
#include "FibonacciGrid.h"

namespace NCXC {

std::vector<Point> fibonacci_sphere(int samples) {
    std::vector<Point> points;
    double phi = M_PI * (std::sqrt(5.0) - 1.0); // golden angle in radians

    for (int i = 0; i < samples; ++i) {
        double y = 1.0 - (i / double(samples - 1)) * 2.0; // y goes from 1 to -1
        double radius = std::sqrt(1.0 - y * y); // radius at y

        double theta = phi * i; // golden angle increment

        double x = std::cos(theta) * radius;
        double z = std::sin(theta) * radius;

        double w = 1.0 / samples; 

        points.push_back({x, y, z, w});
    }

    return points;
}}
/*
int main() {
    int samples = 10;
    auto points = fibonacci_sphere(samples);

    double sum_weights = 0.0;
    for (const auto& p : points) {
        std::cout << "x: " << p.x << ", y: " << p.y 
                  << ", z: " << p.z << ", w: " << p.w << std::endl;
        sum_weights += p.w;
    }

    std::cout << "Sum of weights: " << sum_weights << std::endl;

    return 0;
}
*/