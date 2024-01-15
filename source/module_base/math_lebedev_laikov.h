#ifndef MATH_LEBEDEV_LAIKOV_H
#define MATH_LEBEDEV_LAIKOV_H

#include "vector3.h"
#include <set>
#include <string>

namespace ModuleBase
{

class Lebedev_laikov_grid
{
public:
    Lebedev_laikov_grid(int degree);
    ~Lebedev_laikov_grid();


    void generate_grid_points();

    const ModuleBase::Vector3<double>* get_grid_coor() const
    {
        return grid_coor;
    };

    const double* get_weight() const
    {
        return weight;
    };

    void print_grid_and_weight(std::string filename);

    // degree: can only take the following values
    // degree = { 6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266, 302, 350, 434, 590, 770, 974, 
    // 1202, 1454, 1730, 2030, 2354, 2702, 3074, 3470, 3890, 4334, 4802, 5294, 5810};
    int degree = 6;

private:
    int getLebedevReccurencePoints(int type, int start, double a, double b, double v);

    std::set<int> allowed_degree = { 
        6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266, 
        302, 350, 434, 590, 770, 974, 1202, 1454, 1730, 2030, 
        2354, 2702, 3074, 3470, 3890, 4334, 4802, 5294, 5810
    };

    ModuleBase::Vector3<double> *grid_coor = nullptr;
    double* weight = nullptr;
};

}

#endif