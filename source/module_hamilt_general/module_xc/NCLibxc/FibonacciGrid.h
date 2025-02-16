// FibonacciGrid.h
#ifndef FIBONACCIGRID_H
#define FIBONACCIGRID_H

#include <vector>

struct Point {
    double x, y, z, w;
};

std::vector<Point> fibonacci_sphere(int samples);

#endif // FIBONACCIGRID_H