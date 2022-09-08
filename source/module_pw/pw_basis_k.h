#ifndef PWBASISK_H
#define PWBASISK_H

#include "pw_basis.h"
namespace ModulePW
{


/**
 * @brief Special pw_basis class. It includes different k-points.
 * @author qianrui, Sunliang on 2021-10-15
 * @details
 * Math:
 * plane waves: <r|g,k> = 1/sqrt(V) * exp(i(k+g)r)
 * f(r) = 1/sqrt(V) * \sum_g{c(g)*exp(i(k+g)r)}
 * c(g,k) = \int f(r)*exp(-i(g+k)r) dr
 * 
 * USAGE：
 * ModulePW::PW_Basis_K pwtest;
 * 0. init mpi for PW_Basis
 * pwtest.inimpi(nproc_in_pool,rank_in_pool,POOL_WORLD);
 * 1. setup FFT grids for PW_Basis
 * pwtest.initgrids(lat0,latvec,gridecut);
 * pwtest.initgrids(lat0,latvec,N1,N2,N3); 
 * //double lat0：unit length, (unit: bohr)
 * //ModuleBase::Matrix3 latvec：lattice vector, (unit: lat0), e.g. ModuleBase::Matrix3 latvec(1, 1, 0, 0, 2, 0, 0, 0, 2);
 * //double gridecut：cutoff energy to generate FFT grids, (unit: Ry)
 * //int N1,N2,N3: FFT grids
 * 2. init parameters
 * pwtest.initparameters(gamma_only, ggecut, nks, kvec_d, dividemthd);
 * //bool gamma_only: if use gamma_only
 * //double ggecut: cutoff kinetic energy for planewaves(unit in Ry) (G+K)^2 < ggecut
 * //int nks: number of k points in current cores
 * //ModuleBase::Vector<double>* kvec_d: different k points
 * //int dividemthd: method to divide planewaves to different cores
 * 3. Setup transforms from real space to reciprocal space or from reciprocal space to real space.
 * pwtest.setuptransform(); 
 * pwtest.recip2real(wfg,wfr,ik); //wfg to wfr
 * pwtest.real2recip(wfr,wfg,ik); //wfr to wfg
 * 4. Generate the wave vector for planewaves
 * pwtest.collect_local_pw(); 
 * //then we can use pwtest.gk2, pwtest.gcar, (unit in lat0^-1 or lat0^-2)
 * //getgk2(ik,ig) : get pwtest.gk2: (G+K)^2
 * //getgcar(ik,ig): get pwtest.gcar: G
 * //getgdirect(ik,ig)：get pwtest.gcar: latvec * G
 * //getgpluskcar(ik.ig):   get G+K
 * //getigl2isz(ik,ig): get pwtest.igl2isz_k
 * //getigl2ig(ik,ig):  get pwtest.igl2ig_k
 * 
 */
class PW_Basis_K : public PW_Basis
{

public:
    PW_Basis_K();
    ~PW_Basis_K();

    //init parameters of pw_basis_k class
    void initparameters(
        const bool gamma_only_in,
        const double ecut_in,
        const int nk_in, //number of k points in this pool
        const ModuleBase::Vector3<double> *kvec_d, // Direct coordinates of k points
        const int distribution_type_in = 1,
        const bool xprime_in = true
    );


public:
    int nks=0;//number of k points in this pool
    ModuleBase::Vector3<double> *kvec_d=nullptr; // Direct coordinates of k points
    ModuleBase::Vector3<double> *kvec_c=nullptr; // Cartesian coordinates of k points
    int *npwk=nullptr; //[nks] number of plane waves of different k-points
    int npwk_max=0; //max npwk among all nks k-points, it may be smaller than npw
                  //npw cutoff: (|g|+|k|)^2, npwk in the the npw ball, thus is smaller
    double gk_ecut=0; //Energy cut off for (g+k)^2/2

public:
    //prepare for transforms between real and reciprocal spaces
    void setuptransform();

    int *igl2isz_k=nullptr; //[npwk_max*nks] map (igl,ik) to (is,iz) 
    int *igl2ig_k=nullptr;//[npwk_max*nks] map (igl,ik) to ig

    double *gk2=nullptr; // modulus (G+K)^2 of G vectors [npwk_max*nks]

    //collect gdirect, gcar, gg
    void collect_local_pw();

private:
    //create igl2isz_k map array for fft
    void setupIndGk();
    //calculate G+K, it is a private function
    ModuleBase::Vector3<double> cal_GplusK_cartesian(const int ik, const int ig) const;

public:
    void real2recip(const double* in, std::complex<double>* out, const int ik, const bool add = false, const double factor = 1.0); //in:(nplane,nx*ny)  ; out(nz, ns)
    void real2recip(const std::complex<double>* in, std::complex<double>* out, const int ik, const bool add = false, const double factor = 1.0); //in:(nplane,nx*ny)  ; out(nz, ns)
    void recip2real(const std::complex<double>* in, double* out, const int ik, const bool add = false, const double factor = 1.0); //in:(nz, ns)  ; out(nplane,nx*ny)
    void recip2real(const std::complex<double>* in, std::complex<double> * out, const int ik, const bool add = false, const double factor = 1.0); //in:(nz, ns)  ; out(nplane,nx*ny)

#ifdef __MIX_PRECISION
    void real2recip(const float* in, std::complex<float>* out, const int ik, const bool add = false, const float factor = 1.0); //in:(nplane,nx*ny)  ; out(nz, ns)
    void real2recip(const std::complex<float>* in, std::complex<float>* out, const int ik, const bool add = false, const float factor = 1.0); //in:(nplane,nx*ny)  ; out(nz, ns)
    void recip2real(const std::complex<float>* in, float* out, const int ik, const bool add = false, const float factor = 1.0); //in:(nz, ns)  ; out(nplane,nx*ny)
    void recip2real(const std::complex<float>* in, std::complex<float>* out, const int ik, const bool add = false, const float factor = 1.0); //in:(nz, ns)  ; out(nplane,nx*ny)
#endif

public:
    //operator:
    //get (G+K)^2:
    double& getgk2(const int ik, const int igl) const;
    //get G
    ModuleBase::Vector3<double>& getgcar(const int ik, const int igl) const;
    //get G-direct
    ModuleBase::Vector3<double> getgdirect(const int ik, const int igl) const;
    //get (G+K)
    ModuleBase::Vector3<double> getgpluskcar(const int ik, const int igl) const;
    //get igl2isz_k
    int& getigl2isz(const int ik, const int igl) const;
    //get igl2ig_k or igk(ik,ig) in older ABACUS
    int& getigl2ig(const int ik, const int igl) const;
};

}
#endif //PlaneWave_K class

#include "./pw_basis_k_big.h" //temporary it will be removed

