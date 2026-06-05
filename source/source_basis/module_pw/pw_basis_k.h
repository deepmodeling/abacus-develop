#ifndef PWBASISK_H
#define PWBASISK_H

#include "pw_basis.h"
#include "source_base/module_device/device.h"
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
 * USAGE:
 * ModulePW::PW_Basis_K pwtest;
 * 0. init mpi for PW_Basis
 * pwtest.inimpi(nproc_in_pool,rank_in_pool,POOL_WORLD);
 * 1. setup FFT grids for PW_Basis
 * pwtest.initgrids(lat0,latvec,gridecut);
 * pwtest.initgrids(lat0,latvec,N1,N2,N3);
 * //double lat0: unit length, (unit: bohr)
 * //ModuleBase::Matrix3 latvec: lattice vector, (unit: lat0), e.g. ModuleBase::Matrix3 latvec(1, 1, 0, 0, 2, 0, 0, 0,
 * 2);
 * //double gridecut: cutoff energy to generate FFT grids, (unit: Ry)
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
 * // double erf_ecut_in: the value of the constant energy cutoff
 * // double erf_height_in: the height of the energy step for reciprocal vectors
 * // double erf_sigma_in: the width of the energy step for reciprocal vectors
 * //then we can use pwtest.gk2, pwtest.gcar, (unit in lat0^-1 or lat0^-2)
 * //getgk2(ik,ig) : get pwtest.gk2: (G+K)^2
 * //getgcar(ik,ig): get pwtest.gcar: G
 * //getgdirect(ik,ig): get pwtest.gcar: latvec * G
 * //getgpluskcar(ik.ig):   get G+K
 * //getigl2isz(ik,ig): get pwtest.igl2isz_k
 * //getigl2ig(ik,ig):  get pwtest.igl2ig_k
 *
 */
class PW_Basis_K : public PW_Basis
{

public:
    struct KCacheStats : public PW_Basis::CacheStats
    {
        std::uint64_t gcar_hits = 0;
        std::uint64_t gcar_misses = 0;
        std::uint64_t gk2_hits = 0;
        std::uint64_t gk2_misses = 0;
    };

    PW_Basis_K();
    PW_Basis_K(std::string device_, std::string precision_) : PW_Basis(device_, precision_) {classname="PW_Basis_K";}
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

    int *igl2isz_k=nullptr, * d_igl2isz_k = nullptr; //[npwk_max*nks] map (igl,ik) to (is,iz)
    int *igl2ig_k=nullptr;//[npwk_max*nks] map (igl,ik) to ig
    int *ig2ixyz_k=nullptr; ///< [npw] map ig to ixyz
    std::vector<int> ig2ixyz_k_cpu; /// [npw] map ig to ixyz,which is used in dsp fft.
    double *gk2=nullptr; // modulus (G+K)^2 of G vectors [npwk_max*nks]

    // liuyu add 2023-09-06
    double erf_ecut=0.0;   // the value of the constant energy cutoff
    double erf_height=0.0; // the height of the energy step for reciprocal vectors
    double erf_sigma=0.0;  // the width of the energy step for reciprocal vectors

    //collect gdirect, gcar, gg
    void collect_local_pw(const double& erf_ecut_in = 0.0,
                          const double& erf_height_in = 0.0,
                          const double& erf_sigma_in = 0.1);

    KCacheStats get_k_cache_stats() const;
    void reset_k_cache_stats();

  private:
    void invalidate_cache() override
    {
      PW_Basis::invalidate_cache();
      this->gcar_cache_valid.store(false);
      this->gk_cache_valid.store(false);
      this->gk2 = nullptr;
    }

    void clear_k_cache_storage();
    void sync_gcar_device_cache();
    void sync_gk2_device_cache();

    std::atomic<bool> gcar_cache_valid{false};
    std::atomic<bool> gk_cache_valid{false};
    std::unique_ptr<ModuleBase::Vector3<double>[]> k_gcar_cache_storage;
    std::unique_ptr<double[]> k_gk2_cache_storage;
    std::atomic<std::uint64_t> gcar_cache_hits{0};
    std::atomic<std::uint64_t> gcar_cache_misses{0};
    std::atomic<std::uint64_t> gk2_cache_hits{0};
    std::atomic<std::uint64_t> gk2_cache_misses{0};
    float  * s_gk2 = nullptr;
    double * d_gk2 = nullptr; // modulus (G+K)^2 of G vectors [npwk_max*nks]
    //create igl2isz_k map array for fft
    void setupIndGk();
    // get ig2ixyz_k
    void get_ig2ixyz_k();
  public:
    template <typename FPTYPE>
    void real2recip(const FPTYPE* in,
                    std::complex<FPTYPE>* out,
                    const int ik,
                    const bool add = false,
                    const FPTYPE factor = 1.0) const; // in:(nplane,nx*ny)  ; out(nz, ns)
    template <typename FPTYPE>
    void real2recip(const std::complex<FPTYPE>* in,
                    std::complex<FPTYPE>* out,
                    const int ik,
                    const bool add = false,
                    const FPTYPE factor = 1.0) const; // in:(nplane,nx*ny)  ; out(nz, ns)
    template <typename FPTYPE>
    void recip2real(const std::complex<FPTYPE>* in,
                    FPTYPE* out,
                    const int ik,
                    const bool add = false,
                    const FPTYPE factor = 1.0) const; // in:(nz, ns)  ; out(nplane,nx*ny)
    template <typename FPTYPE>
    void recip2real(const std::complex<FPTYPE>* in,
                    std::complex<FPTYPE>* out,
                    const int ik,
                    const bool add = false,
                    const FPTYPE factor = 1.0) const; // in:(nz, ns)  ; out(nplane,nx*ny)
    #if defined(__DSP)
    template <typename FPTYPE, typename Device>
    void convolution(const Device* ctx,
                      const int ik,
                      const int size,
                      const std::complex<FPTYPE>* input,
                      const FPTYPE*               input1,
                      std::complex<FPTYPE>*       output,
                      const bool add = false,
                      const FPTYPE factor =1.0) const ;

    template <typename FPTYPE>
    void real2recip_dsp(const std::complex<FPTYPE>* in,
                       std::complex<FPTYPE>* out,
                       const int ik,
                       const bool add = false,
                       const FPTYPE factor = 1.0) const; // in:(nplane,nx*ny)  ; out(nz, ns)
    template <typename FPTYPE>
    void recip2real_dsp(const std::complex<FPTYPE>* in,
                       std::complex<FPTYPE>* out,
                       const int ik,
                       const bool add = false,
                       const FPTYPE factor = 1.0) const; // in:(nz, ns)  ; out(nplane,nx*ny)
    
    #endif

     template <typename FPTYPE, typename Device>
    void real_to_recip(const Device* ctx,
                       const std::complex<FPTYPE>* in,
                       std::complex<FPTYPE>* out,
                       const int ik,
                       const bool add = false,
                       const FPTYPE factor = 1.0) const; // in:(nplane,nx*ny)  ; out(nz, ns)
    template <typename FPTYPE, typename Device>
    void recip_to_real(const Device* ctx,
                       const std::complex<FPTYPE>* in,
                       std::complex<FPTYPE>* out,
                       const int ik,
                       const bool add = false,
                       const FPTYPE factor = 1.0) const; // in:(nz, ns)  ; out(nplane,nx*ny)


    template <typename TK,
              typename Device,
              typename std::enable_if<std::is_same<Device, base_device::DEVICE_CPU>::value, int>::type = 0>
    void real_to_recip(const TK* in,
                       TK* out,
                       const int ik,
                       const bool add = false,
                       const typename GetTypeReal<TK>::type factor = 1.0) const
    {
      #if defined(__DSP)
        this->real2recip_dsp(in, out, ik, add, factor);
      #else
        this->real2recip(in,out,ik,add,factor);
      #endif
    }
    template <typename TK,
              typename Device,
              typename std::enable_if<std::is_same<Device, base_device::DEVICE_CPU>::value, int>::type = 0>
    void recip_to_real(const TK* in,
                       TK* out,
                       const int ik,
                       const bool add = false,
                       const typename GetTypeReal<TK>::type factor = 1.0) const
    {
      
      #if defined(__DSP)
        this->recip2real_dsp(in,out,ik,add,factor);
      #else
        this->recip2real(in,out,ik,add,factor);
      #endif
    }
    template <typename FPTYPE>
    void real2recip_gpu(const std::complex<FPTYPE>* in,
                    std::complex<FPTYPE>* out,
                    const int ik,
                    const bool add = false,
                    const FPTYPE factor = 1.0) const; // in:(nplane,nx*ny)  ; out(nz, ns)
                    
    template <typename FPTYPE>
    void recip2real_gpu(const std::complex<FPTYPE>* in,
                    std::complex<FPTYPE>* out,
                    const int ik,
                    const bool add = false,
                    const FPTYPE factor = 1.0) const; // in:(nz, ns)  ; out(nplane,nx*ny)

    template <typename FPTYPE,
              typename Device,
              typename std::enable_if<!std::is_same<Device, base_device::DEVICE_CPU>::value, int>::type = 0>
    void real_to_recip(const FPTYPE* in,
                       FPTYPE* out,
                       const int ik,
                       const bool add = false,
                       const typename GetTypeReal<FPTYPE>::type factor = 1.0) const
    {
        this->real2recip_gpu(in, out, ik, add, factor);
    }

    template <typename TK,
              typename Device,
              typename std::enable_if<std::is_same<Device, base_device::DEVICE_GPU>::value, int>::type = 0>
    void recip_to_real(const TK* in,
                       TK* out,
                       const int ik,
                       const bool add = false,
                       const typename GetTypeReal<TK>::type factor = 1.0) const
    {
        this->recip2real_gpu(in, out, ik, add, factor);
    }

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

    //get ig_to_ix
    std::vector<int> get_ig2ix(const int ik) const;
    //get ig_to_iy
    std::vector<int> get_ig2iy(const int ik) const;
    //get ig_to_iz
    std::vector<int> get_ig2iz(const int ik) const;

    template <typename FPTYPE> FPTYPE * get_gk2_data() const;
    template <typename FPTYPE> FPTYPE * get_gcar_data() const;
    template <typename FPTYPE> FPTYPE * get_kvec_c_data() const;

private:
    float * s_gcar = nullptr, * s_kvec_c = nullptr;
    double * d_gcar = nullptr, * d_kvec_c = nullptr;
};

}
#endif //PlaneWave_K class

#include "./pw_basis_k_big.h" //temporary it will be removed
