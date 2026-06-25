//==========================================================
// AUTHOR : mohan
// LAST UPDATE : 2009-07-07
//==========================================================
#ifndef PW_COMPLEMENT_H
#define PW_COMPLEMENT_H

#include "../src_spillage/tools.h"

namespace PW_complement
{
int get_total_pw_number(const double& ggcut_start, const double& ggcut_end,
                        const int& nx, const int& ny, const int& nz, const Matrix3& GGT);

void get_total_pw(double* ig2sort, Vector3<double> *igsort,
                  const double& ggcut_start, const double& ggcut_end,
                  const int& nx, const int& ny, const int& nz, const Matrix3& GGT, // GGT = G*GT
                  int& ngm); // number of total plane waves.

// output nx, ny, nz according to input: latve and ggcut.
void get_FFT_dimension(const Matrix3& latvec, const double &ggcut, int &nx, int &ny, int &nz);

//==========================================================
// MEMBER FUNCTION :
// NAME : PW_Basis::setup_GVectors
// Second part of the initialization : find out all G
// vectors that |G|^2<=G2max and map it into a one
// dimentional array G1d in the increase order of |G|^2.
// Next generate the indices between the 1D array and
// the 3D G-grid and the FFT grid.
// generate : gg_global, g_global, ig_global
//==========================================================
void setup_GVectors(const Matrix3& G, const int &ngmc_g, double* gg,
                    Vector3<double>* ig, Vector3<double>* g);
}


#endif

