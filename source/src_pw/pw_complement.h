//==========================================================
// AUTHOR : mohan
// LAST UPDATE : 2009-07-07
//==========================================================
#ifndef PW_COMPLEMENT_H
#define PW_COMPLEMENT_H

#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "../module_base/matrix3.h"

namespace PW_complement
{
int get_total_pw_number(const double& ggcut_start, const double& ggcut_end,
                        const int& nx, const int& ny, const int& nz, const ModuleBase::Matrix3& GGT);

void get_total_pw(double* ig2sort, ModuleBase::Vector3<double> *igsort,
                  const double& ggcut_start, const double& ggcut_end,
                  const int& nx, const int& ny, const int& nz, const ModuleBase::Matrix3& GGT, // GGT = G*GT
                  int& ngm); // number of total plane waves.

// output nx, ny, nz according to input: latve and ggcut.
void get_FFT_dimension(const ModuleBase::Matrix3& latvec, const double &ggcut, int &nx, int &ny, int &nz,
const int &bx, const int &by, const int &bz);

//==========================================================
// MEMBER FUNCTION :
// NAME : PW_Basis::setup_GVectors
// Second part of the initialization : find out all G
// vectors that |G|^2<=G2max and std::map it into a one
// dimentional array G1d in the increase order of |G|^2.
// Next generate the indices between the 1D array and
// the 3D G-grid and the FFT grid.
// generate : gg_global, g_global, ig_global
//==========================================================
void setup_GVectors(const ModuleBase::Matrix3& G, const int &ngmc_g, double* gg,
                    ModuleBase::Vector3<double>* ig, ModuleBase::Vector3<double>* g);

void get_total_pw_after_vc(double* ig2sort0, double* ig2sort, ModuleBase::Vector3<double> *igsort,
                  const double& ggcut_start, const double& ggcut_end,
                  const int& nx, const int& ny, const int& nz, const ModuleBase::Matrix3& GGT, // GGT = G*GT
                  const ModuleBase::Matrix3& GGT0,
                  int& ngm); // number of total plane waves.

// #ifndef __MPI

void get_ngmw(const int &ngmc, const double& ggwfc2, const double* gg_global, int &ngmw);

void get_ig2fftw(const int &ngmw, const int &nx, const int &ny, const int &nz,
                 const ModuleBase::Vector3<double> *ig, int *ig2fftw);

//==========================================================
// (1) allocate ig2fftc, ig1, ig2, ig3
// (2) use ig and ncx, ncy, ncz to give the arrays value.
//==========================================================
void get_ig2fftc(const int &ngmc, const int &ncx, const int &ncy, const int &ncz,
                 const ModuleBase::Vector3<double> *ig, int* ig1, int *ig2, int *ig3, int *ig2fftc);
// #endif
}


#endif

