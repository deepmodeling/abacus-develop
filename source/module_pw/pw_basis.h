#ifndef PWBASIS_H
#define PWBASIS_H

#include "../module_base/matrix.h"
#include "../module_base/matrix3.h"
#include "../module_base/vector3.h"
#include <complex>
#include "fft.h"

namespace ModulePW
{

/**
 * @brief A class which can convert a function of "r" to the corresponding linear 
 * superposition of plane waves (real space to reciprocal space)
 * or convert a linear superposition of plane waves to the function 
 * of "r" (reciprocal to real).
 * @author qianrui, Sunliang on 2021-10-15
 * @details
 * Math:
 * plane waves: <r|g>=1/sqrt(V) * exp(igr)
 * f(r) = 1/sqrt(V) * \sum_g{c(g)*exp(igr)}
 * c(g) = \int f(r)*exp(-igr) dr
 * USAGE：
 * ModulePW::PW_Basis pwtest;
 * 1. setup FFT grids for PW_Basis
 * pwtest.initgrids(lat0,latvec,gridecut,nproc_in_pool,rank_in_pool);
 * pwtest.initgrids(lat0,latvec,N1,N2,N3,nproc_in_pool,rank_in_pool); 
 * //double lat0：unit length, (unit: bohr)
 * //ModuleBase::Matrix3 latvec：lattice vector, (unit: lat0), e.g. ModuleBase::Matrix3 latvec(1, 1, 0, 0, 2, 0, 0, 0, 2);
 * //double gridecut：cutoff energy to generate FFT grids, (unit: Ry)
 * //int N1,N2,N3: FFT grids
 * 2. init parameters
 * pwtest.initparameters(gamma_only,ggecut,dividemthd);
 * //bool gamma_only: if use gamma_only
 * //double ggecut: cutoff kinetic energy for planewaves,(unit in Ry) G^2 < ggecut
 * //int dividemthd: method to divide planewaves to different cores
 * 3. Setup transforms from real space to reciprocal space or from reciprocal space to real space.
 * pwtest.setuptransform(); 
 * pwtest.recip2real(rhog,rhor); //rhog to rhor
 * pwtest.real2recip(rhor,rhog); //rhor to rhog
 * 4. Generate the wave vector for planewaves
 * pwtest.collect_local_pw(); 
 * //then we can use pwtest.gg, pwtest.gdirect, pwtest.gcar, (unit in lat0^-1 or lat0^-2)
 * 
 */
class PW_Basis
{

public:
    PW_Basis();
    virtual ~PW_Basis();

    //Init the grids for FFT
    virtual void initgrids(
        const double lat0_in, //unit length (unit in bohr)
        const ModuleBase::Matrix3 latvec_in, // Unitcell lattice vectors (unit in lat0) 
        const double gridecut, //unit in Ry, ecut to set up grids
        const int poolnproc_in, // Number of processors in this pool
        const int poolrank_in // Rank in this pool
    );
    //Init the grids for FFT
    virtual void initgrids(
        const double lat0_in,
        const ModuleBase::Matrix3 latvec_in, // Unitcell lattice vectors
        const int nx_in, int ny_in, int nz_in,
        const int poolnproc_in, // Number of processors in this pool
        const int poolrank_in // Rank in this pool
    );

    //Init some parameters
    void initparameters(
        const bool gamma_only_in,
        const double pwecut_in, //unit in Ry, ecut to decides plane waves
        const int distribution_type_in = 1
    );

//===============================================
//                 distribution maps
//===============================================
public:
    //reciprocal-space
    // only on first proc.
    
    int *ig2isz=nullptr; // map ig to (is, iz).
    int *istot2ixy=nullptr; // istot2ixy[is]: iy + ix * ny of is^th stick among all sticks.
    int *is2fftixy=nullptr; // is2fftixy[is]: iy + ix * ny of is^th stick among sticks on current proc.
    int *fftixy2ip=nullptr; // fftixy2ip[iy + ix * fftny]: ip of proc which contains stick on (ix, iy). if no stick: -1
    int nst=0; //num. of sticks in current proc.
    int *nst_per=nullptr;// nst on each core
    int nstnz=0; // nst * nz
    int nstot=0; //num. of sticks in total.
    int npw=0; //num. of plane waves in current proc.
    int *npw_per=nullptr; //npw on each core
    int npwtot=0; // total num. of plane waves in all proc. in this pool

    //real space
    int nrxx=0; //num. of real space grids
    int *startz=nullptr; //startz[ip]: starting z plane in the ip-th proc. in current POOL_WORLD 
	int *numz=nullptr; //numz[ip]: num. of z planes in the ip-th proc. in current POOL_WORLD
    int *numg=nullptr; //numg[ip] :  nst_per[poolrank] * numz[ip] 
    int *numr=nullptr; //numr[ip] :  numz[poolrank] * nst_per[ip]
    int *startg=nullptr;  // startg[ip] = numg[ip-1] + startg[ip-1]
    int *startr=nullptr;  // startr[ip] = numr[ip-1] + startr[ip-1]
    int startz_current=0;
    int nplane=0; //num. of planes in current proc.

    ModuleBase::Vector3<double> *gdirect=nullptr;		//(= *G1d) ; // ig = new Vector igc[npw]
    ModuleBase::Vector3<double> *gcar=nullptr;   			//G vectors in cartesian corrdinate
    double *gg=nullptr;       	// modulus (G^2) of G vectors [npw]
    //gg[ng]=ig[ng]*GGT*ig[ng]/(lat0*lat0)=g[ng]*g[ng] (/lat0*lat0)
	// gg_global dimension: [cutgg_num_now] (save memory skill is used)
    int ig_gge0=-1;    //ig when gg == 0

    //distribute plane waves and grids and set up fft
    void setuptransform();

protected:
    int *startnsz_per=nullptr;//useless intermediate variable// startnsz_per[ip]: starting is * nz stick in the ip^th proc.

    //distribute plane waves to different processors
    void distribute_g();

    //distribute real-space grids to different processors
    virtual void distribute_r();

    //prepare for MPI_Alltoall
    void getstartgr();


public:
    //prepare for transforms between real and reciprocal spaces
    void collect_local_pw();

    // void collect_tot_pw(
    //     double* gg_global,
    //     ModuleBase::Vector3<double> *gdirect_global,
    //     ModuleBase::Vector3<double> *gcar_global
    // ); 
public:
    int ngg=0; //number of different modulus (G^2) of G vectors
    int *ig2igg=nullptr;//[npw] map ig to igg(<ngg: the index of G^2)
    double *gg_uniq=nullptr; //[ngg] modulus (G^2) of G vectors of igg, each gg of igg is unique.
    void collect_uniqgg();
   

public:
    bool gamma_only=false;	// only half g are used.
    double ggecut=0;    //Energy cut off for g^2/2, unit in 1/lat0^2, ggecut=ecutwfc(Ry)*lat0^2/4pi^2
    double lat0=1;     //unit length for lattice, unit in bohr
    double tpiba=1;    //  2pi/lat0
    double tpiba2=1;   //  4pi^2/lat0^2
    ModuleBase::Matrix3 latvec; // Unitcell lattice vectors, unit in lat0
    ModuleBase::Matrix3 G; // reciprocal lattice vector, unit in 1/lat0
    ModuleBase::Matrix3 GT; // traspose of G
    ModuleBase::Matrix3 GGT; // GGT = G*GT
    int distribution_type=1;
    int poolnproc=1;
    int poolrank=0;
    //distribute plane waves to different processors

protected:
    //method 1: first consider number of plane waves
    void distribution_method1(); 
    // Distribute sticks to cores in method 1.
    void divide_sticks_1(
        int* st_i,          // x or x + fftnx (if x < 0) of stick.
        int* st_j,          // y or y + fftny (if y < 0) of stick.
        int* st_length     // the stick on (ix, iy) consists of st_length[ix*fftny+iy] planewaves.
    );

    //method 2: first consider number of sticks
    void distribution_method2();
    // Distribute sticks to cores in method 2.
    void divide_sticks_2();
   
    //Count the total number of planewaves (tot_npw) and sticks (this->nstot) (in distributeg method1 and method2)
    void count_pw_st(
        int* st_length2D, // the number of planewaves that belong to the stick located on (x, y).
        int* st_bottom2D  // the z-coordinate of the bottom of stick on (x, y).
    );

    //get ig2isz and is2fftixy
    void get_ig2isz_is2fftixy(
        int* st_bottom,     // minimum z of stick, stored in 1d array with tot_nst elements.
        int* st_length     // the stick on (x, y) consists of st_length[x*fftny+y] planewaves.
    );

    //Collect the x, y indexs, length of the sticks (in distributeg method1)
    void collect_st(
        int* st_length2D,                               // the number of planewaves that belong to the stick located on (x, y), stored in 2d x-y plane.
        int* st_bottom2D,                               // the z-coordinate of the bottom of stick on (x, y), stored in 2d x-y plane.
        int* st_i,                                      // x or x + fftnx (if x < 0) of stick.
        int* st_j,                                      // y or y + fftny (if y < 0) of stick.
        int* st_length                                 // number of planewaves in stick, stored in 1d array with tot_nst elements.
    );

    //get istot2ixy
    void get_istot2ixy(
        int* st_i,          // x or x + fftnx (if x < 0) of stick.
        int* st_j          // y or y + fftny (if y < 0) of stick.
    );

    //Create the maps from ixy to (in method 2)
    void create_maps(
        int* st_length2D  // the number of planewaves that belong to the stick located on (x, y), stored in 2d x-y plane.
    );

//===============================================
//                  FFT
//===============================================
public:
	// FFT dimensions for wave functions.
	int fftnx=0, fftny=0, fftnz=0, fftnxyz=0, fftnxy=0;
    int nx=0, ny=0, nz=0, nxyz=0, nxy=0; // Gamma_only: fftny = int(ny/2)-1 , others: fftny = ny
    int liy=0,riy=0;// liy: the left edge of the pw ball; riy: the right edge of the pw ball
    int nmaxgr=0; // Gamma_only: max between npw and (nrxx+1)/2, others: max between npw and nrxx
                // Thus complex<double>[nmaxgr] is able to contain either reciprocal or real data
    FFT ft;
    //The position of pointer in and out can be equal(in-place transform) or different(out-of-place transform).
    void real2recip(const double * in, std::complex<double> * out, const bool add = false, const double factor = 1.0); //in:(nplane,nx*ny)  ; out(nz, ns)
    void real2recip(const std::complex<double> * in, std::complex<double> * out, const bool add = false, const double factor = 1.0); //in:(nplane,nx*ny)  ; out(nz, ns)
    void recip2real(const std::complex<double> * in, double *out, const bool add = false, const double factor = 1.0); //in:(nz, ns)  ; out(nplane,nx*ny)
    void recip2real(const std::complex<double> * in, std::complex<double> * out, const bool add = false, const double factor = 1.0); //in:(nz, ns)  ; out(nplane,nx*ny)

#ifdef __MIX_PRECISION
    void real2recip(const float * in, std::complex<float> * out, const bool add = false, const float factor = 1.0); //in:(nplane,nx*ny)  ; out(nz, ns)
    void real2recip(const std::complex<float> * in, std::complex<float> * out, const bool add = false, const float factor = 1.0); //in:(nplane,nx*ny)  ; out(nz, ns)
    void recip2real(const std::complex<float> * in, float *out, const bool add = false, const float factor = 1.0); //in:(nz, ns)  ; out(nplane,nx*ny)
    void recip2real(const std::complex<float> * in, std::complex<float> * out, const bool add = false, const float factor = 1.0); //in:(nz, ns)  ; out(nplane,nx*ny)
#endif
protected:
    //gather planes and scatter sticks of all processors
    template<typename T>
    void gatherp_scatters(std::complex<T> *in, std::complex<T> *out); 

    //gather sticks of and scatter planes of all processors
    template<typename T>
    void gathers_scatterp(std::complex<T> *in, std::complex<T> *out); 
public:
    //get fftixy2is;
    void getfftixy2is(int * fftixy2is);
};

}
#endif //PlaneWave 

#include "./pw_basis_big.h" //temporary it will be removed