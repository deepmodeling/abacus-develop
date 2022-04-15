#ifndef PWBASIS_H
#define PWBASIS_H

#include "../module_base/matrix.h"
#include "../module_base/matrix3.h"
#include "../module_base/vector3.h"
#include <complex>
#include "fft.h"

namespace ModulePW
{

//
//A class which can convert a function of "r" to the corresponding linear 
// superposition of plane waves (real space to reciprocal space)
// or convert a linear superposition of plane waves to the function 
// of "r" (reciprocal to real).
//plane waves: <r|g>=1/sqrt(V) * exp(igr)
// f(r) = 1/sqrt(V) * \sum_g{c(g)*exp(igr)}
//
class PW_Basis
{

public:
    PW_Basis();
    ~PW_Basis();

    //Init the grids for FFT
    void initgrids(
        double lat0_in, //unit length (unit in bohr)
        ModuleBase::Matrix3 latvec_in, // Unitcell lattice vectors (unit in lat0) 
        double gridecut //unit in Ry, ecut to set up grids
    );
    //Init the grids for FFT
    void initgrids(
        double lat0_in,
        ModuleBase::Matrix3 latvec_in, // Unitcell lattice vectors
        int nx_in, int bigny_in, int nz_in
    );

    //Init some parameters
    void initparameters(
        bool gamma_only_in,
        double pwecut_in, //unit in Ry, ecut to decides plane waves
        int poolnproc_in, // Number of processors in this pool
        int poolrank_in, // Rank in this pool
        int distribution_type_in
    );

//===============================================
//                 distribution maps
//===============================================
public:
    //reciprocal-space
    // only on first proc.
    int *startnsz_per; // startnsz_per[ip]: starting is * nz stick in the ip^th proc.
    int *nstnz_per; // nz * nst(number of sticks) on each core.
    int *nst_per;// nst on each core
    // on all proc.
    int *ig2isz; // map ig to (is, iz).
    int *istot2bigixy; // istot2bigixy[is]: iy + ix * bigny of is^th stick among all sticks.
    int *ixy2istot; // ixy2istot[ix + iy * nx]: is of stick on (ix, iy) among all sticks.
    int *is2ixy; // is2ixy[is]: ix + iy * bignx of is^th stick among sticks on current proc.
    int *ixy2ip; // ixy2ip[ix + iy * nx]: ip of proc which contains stick on (ix, iy).
    int nst; //num. of sticks in current proc.
    int nstnz; // nst * nz
    int nstot; //num. of sticks in total.
    int npw; //num. of plane waves in current proc.
    //real space
    int nrxx; //num. of real space grids
    int *startz; //startz[ip]: starting z plane in the ip-th proc. in current POOL_WORLD 
	int *numz; //numz[ip]: num. of z planes in the ip-th proc. in current POOL_WORLD
    int *numg; //numg[ip] :  nst_per[poolrank] * numz[ip] 
    int *numr; //numr[ip] :  numz[poolrank] * nst_per[ip]
    int *startg;  // startg[ip] = numg[ip-1] + startg[ip-1]
    int *startr;  // startr[ip] = numr[ip-1] + startr[ip-1]
    int nplane; //num. of planes in current proc.

    ModuleBase::Vector3<double> *gdirect;		//(= *G1d) ; // ig = new Vector igc[ngmc]
    ModuleBase::Vector3<double> *gcar;   			//G vectors in cartesian corrdinate
    double *gg;       	// modulus (G^2) of G vectors [ngmc]
    //gg[ng]=ig[ng]*GGT*ig[ng]/(lat0*lat0)=g[ng]*g[ng] (/lat0*lat0)
	// gg_global dimension: [cutgg_num_now] (save memory skill is used)

    //distribute plane waves and grids and set up fft
    void setuptransform();
    
    //distribute plane waves to different processors
    void distribute_g();

    //distribute real-space grids to different processors
    void distribute_r();

    void getstartgr();

    //distribute plane waves to different processors
    void distribution_method1(); // x varies fast
    void distribution_method2(); // sticks sorted according to ixy
    // void distribution_method3(); // y varies fast

    void collect_local_pw();

    // void collect_tot_pw(
    //     double* gg_global,
    //     ModuleBase::Vector3<double> *gdirect_global,
    //     ModuleBase::Vector3<double> *gcar_global
    // ); 
   

public:
    bool gamma_only;	// only half g are used.
    double ggecut;    //Energy cut off for g^2/2, unit in 1/lat0^2, ggecut=ecutwfc(Ry)*lat0^2/4pi^2
    double lat0;     //unit length for lattice, unit in bohr
    ModuleBase::Matrix3 latvec; // Unitcell lattice vectors, unit in lat0
    ModuleBase::Matrix3 G; // reciprocal lattice vector, unit in 1/lat0
    ModuleBase::Matrix3 GT; // traspose of G
    ModuleBase::Matrix3 GGT; // GGT = G*GT
    int distribution_type;
    int poolnproc;
    int poolrank;
   

// for both distributeg_method1 and distributeg_method2
    void count_pw_st(
        int &tot_npw,     // total number of planewaves.
        int* st_length2D, // the number of planewaves that belong to the stick located on (x, y).
        int* st_bottom2D  // the z-coordinate of the bottom of stick on (x, y).
    );
    void get_ig2isz_is2ixy(
        int* st_bottom,     // minimum z of stick, stored in 1d array with tot_nst elements.
        int* st_length     // the stick on (x, y) consists of st_length[x*ny+y] planewaves.
    );
// for distributeg_method1
    void collect_st(
        int* st_length2D,                               // the number of planewaves that belong to the stick located on (x, y), stored in 2d x-y plane.
        int* st_bottom2D,                               // the z-coordinate of the bottom of stick on (x, y), stored in 2d x-y plane.
        int* st_i,                                      // x or x + nx (if x < 0) of stick.
        int* st_j,                                      // y or y + ny (if y < 0) of stick.
        int* st_length                                 // number of planewaves in stick, stored in 1d array with tot_nst elements.
    );
    void divide_sticks(
        int* st_i,          // x or x + nx (if x < 0) of stick.
        int* st_j,          // y or y + ny (if y < 0) of stick.
        int* st_length,     // the stick on (x, y) consists of st_length[x*ny+y] planewaves.
        int* npw_per       // number of planewaves on each core.
    );
    void get_istot2bigixy(
        int* st_i,          // x or x + nx (if x < 0) of stick.
        int* st_j          // y or y + ny (if y < 0) of stick.
    );
// for distributeg_method2
    void divide_sticks2();
    void create_maps(
        int* st_length2D,  // the number of planewaves that belong to the stick located on (x, y), stored in 2d x-y plane.
        int* npw_per       // number of planewaves on each core.
    );
// for distributeg_method3
    // void divide_sticks2(
    //     const int tot_npw,  // total number of planewaves.
    //     int* st_i,          // x or x + nx (if x < 0) of stick.
    //     int* st_j,          // y or y + ny (if y < 0) of stick.
    //     int* st_length,     // the stick on (x, y) consists of st_length[x*ny+y] planewaves.
    //     int* npw_per,       // number of planewaves on each core.
    //     int* nst_per,       // number of sticks on each core.
    //     int* is2ip         // ip of core containing is^th stick, map is to ip.         
    // );
    // void get_istot2ixy2(
    //     int* st_i,          // x or x + nx (if x < 0) of stick.
    //     int* st_j,          // y or y + ny (if y < 0) of stick.
    //     int* is2ip          // ip of core containing is^th stick, map is to ip.
    // );
    // void get_ig2isz_is2ixy2(
    //     int* st_i,          // x or x + nx (if x < 0) of stick.
    //     int* st_j,          // y or y + ny (if y < 0) of stick.
    //     int* st_bottom,     // minimum z of stick, stored in 1d array with tot_nst elements.
    //     int* st_length,     // the stick on (x, y) consists of st_length[x*ny+y] planewaves.
    //     int* is2ip          // ip of core containing is^th stick, map is to ip.
    // );

//===============================================
//                  FFT
//===============================================
public:
	// FFT dimensions for wave functions.
	int nx, ny, nz, nxyz, nxy;
    int bigny, bignxyz, bignxy; // Gamma_only: ny = int(bigny/2)-1 , others: ny = bigny
    int liy,riy;// liy: the left edge of the pw ball; riy: the right edge of the pw ball
    int maxgrids; // max between nz * ns and bignxy * nplane
    FFT ft;

    void real2recip(double * in, std::complex<double> * out); //in:(nplane,nx*ny)  ; out(nz, ns)
    void real2recip(std::complex<double> * in, std::complex<double> * out); //in:(nplane,nx*ny)  ; out(nz, ns)
    void recip2real(std::complex<double> * in, double *out); //in:(nz, ns)  ; out(nplane,nx*ny)
    void recip2real(std::complex<double> * in, std::complex<double> * out); //in:(nz, ns)  ; out(nplane,nx*ny)

#ifdef __MIX_PRECISION
    void real2recip(float * in, std::complex<float> * out); //in:(nplane,nx*ny)  ; out(nz, ns)
    void real2recip(std::complex<float> * in, std::complex<float> * out); //in:(nplane,nx*ny)  ; out(nz, ns)
    void recip2real(std::complex<float> * in, float *out); //in:(nz, ns)  ; out(nplane,nx*ny)
    void recip2real(std::complex<float> * in, std::complex<float> * out); //in:(nz, ns)  ; out(nplane,nx*ny)
#endif
    template<typename T>
    void gatherp_scatters(std::complex<T> *in, std::complex<T> *out); //gather planes and scatter sticks of all processors
    template<typename T>
    void gathers_scatterp(std::complex<T> *in, std::complex<T> *out); //gather sticks of and scatter planes of all processors
    // void gathers_scatterp2(std::complex<double> *in, std::complex<double> *out); //gather sticks of and scatter planes of all processors
    // void gatherp_scatters2(std::complex<double> *in, std::complex<double> *out); //gather sticks of and scatter planes of all processors
    // void gatherp_scatters_gamma(std::complex<double> *in, std::complex<double> *out); //gather planes and scatter sticks of all processors, used when gamma_only
    // void gathers_scatterp_gamma(std::complex<double> *in, std::complex<double> *out); //gather sticks of and scatter planes of all processors, used when gamma only

};

}
#endif //PlaneWave class
