#ifndef MATH_POLYINT_H
#define MATH_POLYINT_H

#include "realarray.h"

using namespace std;

// mohan add 2021-05-07
class PolyInt
{

	public:

	PolyInt();
	~PolyInt();

    //========================================================
    // Polynomial_Interpolation
    //========================================================
    static void Polynomial_Interpolation
    (
        const realArray &table,
        const int &dim1,
        const int &dim2,
        realArray &y,
        const int &dim_y,
        const int &table_length,
        const double &table_interval,
        const double &x
    );

    static double Polynomial_Interpolation
    (
        const realArray &table,
        const int &dim1,
        const int &dim2,
        const int &table_length,
        const double &table_interval,
        const double &x             // input value
    );

    static double Polynomial_Interpolation             // pengfei Li 2018-3-23
    (
        const realArray &table,
        const int &dim1,
        const int &dim2,
        const int &dim3,
        const int &table_length,
        const double &table_interval,
        const double &x             // input value
    );

	static double Polynomial_Interpolation
	(
        const double *table,
        const int &table_length,
        const double &table_interval,
        const double &x             // input value
    );

    static double Polynomial_Interpolation_xy
    (
        const double *xpoint,
        const double *ypoint,
        const int table_length,
        const double &x             // input value
    );

};
#endif
